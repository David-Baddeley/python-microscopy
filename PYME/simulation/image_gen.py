#!/usr/bin/python

##################
# image_gen.py
#
# Cleaner re-implementation of PYME.Acquire.Hardware.Simulator.rend_im, with
# state encapsulated in an ImageGenerator instance rather than module globals,
# and illumination functions as __call__-able classes owning their own state.
#
# rend_im.py is left untouched; this module is a parallel implementation to be
# swapped in incrementally.
#
# Copyright David Baddeley, 2009-2026
# d.baddeley@auckland.ac.nz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##################

import multiprocessing
import threading

import numpy as np
from numpy.fft import fftn, ifftn, ifftshift
from scipy import ndimage

from PYME.Analysis import MetaData
from PYME.Deconv.wiener import resizePSF
from PYME.localization import cInterp


class ConstIllumFunction(object):
    '''Uniform illumination, independent of fluorophore position.'''

    def __call__(self, fluors, position):
        return np.float32(1.0)


class PSFIllumFunction(object):
    '''Illuminate with the (spatially varying) intensity of a PSF model.'''

    def __init__(self, psf_provider, chan=0):
        # psf_provider() -> (model, dx, dy, dz), so a different PSF than the
        # detection PSF can be substituted if desired.
        self.psf_provider = psf_provider
        self.chan = chan

    def __call__(self, fluors, position):
        im, dx, dy, dz = self.psf_provider(self.chan)

        xi = np.maximum(np.minimum(np.round_((fluors['x'] - position[0]) / dx + im.shape[0] / 2).astype('i'), im.shape[0] - 1), 0)
        yi = np.maximum(np.minimum(np.round_((fluors['y'] - position[1]) / dy + im.shape[1] / 2).astype('i'), im.shape[1] - 1), 0)
        zi = np.maximum(np.minimum(np.round_((fluors['z'] - position[2]) / dz + im.shape[2] / 2).astype('i'), im.shape[2] - 1), 0)

        return im[xi, yi, zi]


class ROIIllumFunction(object):
    '''
    Very crude ROI-based illumination. Assumes hard edges, no diffraction.
    '''

    def __init__(self, get_pixelsize_nm, roi_size=256):
        # get_pixelsize_nm() -> (vx, vy) in nm
        self.get_pixelsize_nm = get_pixelsize_nm
        self.roi_size = roi_size

    def __call__(self, fluors, position):
        vx, vy = self.get_pixelsize_nm()

        xi = np.round_((fluors['x'] - position[0]) / vx)
        yi = np.round_((fluors['y'] - position[1]) / vy)

        # cast to f4 - ilFrac is combined with f4 arrays in the illuminate() extension
        return ((xi > 0) * (xi < self.roi_size) * (yi > 0) * (yi < self.roi_size)).astype('f4')


class PatternIllumFunction(object):
    '''FFT-convolved illumination pattern, e.g. for a projected/scanned illumination profile.'''

    def __init__(self, psf_provider):
        # psf_provider() -> (model, dx, dy, dz)
        self.psf_provider = psf_provider

        self.pattern = None
        self.z_offset = 0
        self.dx = self.dy = self.dz = None

        self._cache = None
        self._cache_key = None

    def set_pattern(self, pattern, z0, chan=0):
        '''Convolve `pattern` with the PSF to give a diffraction-limited illumination pattern.'''
        sx, sy = pattern.shape
        im, dx, dy, dz = self.psf_provider(chan)
        psx, psy, sz = im.shape

        il = np.zeros([sx, sy, sz], 'f')
        il[:, :, sz // 2] = pattern
        ps = np.zeros_like(il)
        if sx > psx:
            ps[(sx // 2 - psx // 2):(sx // 2 + psx // 2), (sy // 2 - psy // 2):(sy // 2 + psy // 2), :] = im
        else:
            ps[:, :, :] = im[(psx // 2 - sx // 2):(psx // 2 + sx // 2), (psy // 2 - sy // 2):(psy // 2 + sy // 2), :]
        ps = ps / ps[:, :, sz // 2].sum()

        self.pattern = abs(ifftshift(ifftn(fftn(il) * fftn(ps)))).astype('f')
        self.z_offset = z0
        self.dx, self.dy, self.dz = dx, dy, dz

        self._cache = None
        self._cache_key = None

    def __call__(self, fluors, position):
        key = hash((fluors[0]['x'], fluors[0]['y'], fluors[0]['z']))

        if self._cache is not None and self._cache_key == key:
            return self._cache

        x = fluors['x'] / self.dx + self.pattern.shape[0] / 2
        y = fluors['y'] / self.dy + self.pattern.shape[1] / 2
        z = (fluors['z'] - self.z_offset) / self.dz + self.pattern.shape[2] / 2

        self._cache_key = key
        self._cache = ndimage.map_coordinates(self.pattern, [x, y, z], order=1, mode='nearest')
        return self._cache


class SIMIllumFunction(object):
    '''Sinusoidal structured illumination pattern.'''

    def __init__(self, k=np.pi / 180., theta=0, phi=0):
        self.k = k
        self.theta = theta
        self.phi = phi

    def __call__(self, fluors, position):
        x = fluors['x']
        y = fluors['y']

        kx = np.cos(self.theta) * self.k
        ky = np.sin(self.theta) * self.k

        # cast to f4 - ilFrac is combined with f4 arrays in the illuminate() extension
        return ((1 + np.cos(x * kx + y * ky + self.phi)) / 2).astype('f4')


class ImageGenerator(object):
    '''
    Encapsulates the PSF model, illumination functions, and rendering logic
    needed to simulate camera images from a collection of fluorophores.

    This is a cleaner re-implementation of
    PYME.Acquire.Hardware.Simulator.rend_im, with (former) module-level global
    state as instance attributes.
    '''

    def __init__(self):
        self.IntXVals = None
        self.IntYVals = None
        self.IntZVals = None

        self.interpModel_by_chan = [None, None, None, None]

        self.dx = None
        self.dy = None
        self.dz = None

        self.mdh = MetaData.NestedClassMDHandler(MetaData.TIRFDefault)
        # 50nm z spacing so we are not relying quite as much on interpolation for the z shape
        self.mdh['voxelsize.z'] = 0.05

        self.illumination_functions = {
            'ConstIllum': ConstIllumFunction(),
            'PSFIllumFunction': PSFIllumFunction(self.get_psf_model),
            'ROIIllumFunction': ROIIllumFunction(self._get_pixelsize_nm),
            'patternIllumFcn': PatternIllumFunction(self.get_psf_model),
            'SIMIllumFcn': SIMIllumFunction(),
        }

    def _get_pixelsize_nm(self):
        vs = self.mdh.voxelsize_nm
        return vs.x, vs.y

    def get_psf_model(self, chan=0):
        '''PSF provider contract: returns (model, dx, dy, dz) for the given channel.'''
        im = self.interpModel_by_chan[chan]
        if im is None and chan != 0:
            im = self.interpModel_by_chan[0]

        return im, self.dx, self.dy, self.dz

    def set_pixelsize_nm(self, pixelsize):
        self.mdh['voxelsize.x'] = 1e-3 * pixelsize
        self.mdh['voxelsize.y'] = 1e-3 * pixelsize

    def genTheoreticalModel(self, zernikes={}, **kwargs):
        from PYME.Analysis.PSFGen import fourierHNA

        vs = self.mdh.voxelsize_nm
        self.IntXVals = vs.x * np.mgrid[-150:150]
        self.IntYVals = vs.y * np.mgrid[-150:150]
        self.IntZVals = vs.z * np.mgrid[-30:30]

        self.dx, self.dy, self.dz = vs

        im = fourierHNA.GenZernikeDPSF(self.IntZVals, zernikes, X=self.IntXVals, Y=self.IntYVals, dx=vs.x, **kwargs)

        for i in range(1, len(self.interpModel_by_chan)):
            self.interpModel_by_chan[i] = None

        # normalise to 1 and clip
        self.interpModel_by_chan[0] = np.maximum(im / im[:, :, int(len(self.IntZVals) / 2)].sum(), 0).astype('f4')

    def genTheoreticalModel4Pi(self, zernikes=[{}, {}], phases=[0, np.pi / 2, np.pi, 3 * np.pi / 2], **kwargs):
        from PYME.Analysis.PSFGen import fourierHNA

        vs = self.mdh.voxelsize_nm
        self.IntXVals = vs.x * np.mgrid[-150:150]
        self.IntYVals = vs.y * np.mgrid[-150:150]
        self.IntZVals = 20 * np.mgrid[-60:60]

        self.dx, self.dy = vs.x, vs.y
        self.dz = 20.

        for i, phase in enumerate(phases):
            im = fourierHNA.Gen4PiPSF(self.IntZVals, phi=phase, zernikeCoeffs=zernikes, X=self.IntXVals, Y=self.IntYVals, dx=vs.x, **kwargs)

            zm = int(len(self.IntZVals) / 2)
            # due to interference we can have slices with really low sum
            norm = im[:, :, (zm - 10):(zm + 10)].sum(1).sum(0).max()
            self.interpModel_by_chan[i] = np.maximum(im / norm, 0).astype('f4')

    def get_psf_image_stack(self):
        from PYME.IO.image import ImageStack
        from PYME.IO.MetaDataHandler import NestedClassMDHandler

        mdh = NestedClassMDHandler()
        mdh['ImageType'] = 'PSF'
        mdh['voxelsize.x'] = self.dx / 1e3
        mdh['voxelsize.y'] = self.dy / 1e3
        mdh['voxelsize.z'] = self.dz / 1e3

        return ImageStack(data=[c for c in self.interpModel_by_chan if c is not None], mdh=mdh, titleStub='Simulated PSF')

    def setModel(self, modName):
        from PYME.IO import load_psf

        mod, vs_nm = load_psf.load_psf(modName)
        mod = resizePSF(mod, self.get_psf_model()[0].shape)

        self.IntXVals = vs_nm.x * np.mgrid[-(mod.shape[0] / 2.):(mod.shape[0] / 2.)]
        self.IntYVals = vs_nm.y * np.mgrid[-(mod.shape[1] / 2.):(mod.shape[1] / 2.)]
        self.IntZVals = vs_nm.z * np.mgrid[-(mod.shape[2] / 2.):(mod.shape[2] / 2.)]

        self.dx, self.dy, self.dz = vs_nm

        # normalise to 1 and clip
        self.interpModel_by_chan[0] = np.maximum(mod / mod[:, :, len(self.IntZVals) // 2].sum(), 0).astype('f4')

    def set_illum_pattern(self, pattern, z0):
        self.illumination_functions['patternIllumFcn'].set_pattern(pattern, z0)

    def _render_fluor_subset(self, im, fl, A, x0, y0, z, dx, dy, dz, maxz, ChanXOffsets=[0, ], ChanZOffsets=[0, ], ChanSpecs=None):
        if ChanSpecs is None:
            z_ = np.clip(z - fl['z'], -maxz, maxz).astype('f')
            roiSize = np.minimum(8 + np.abs(z_) * (2.5 / dx), 140).astype('i')
            cInterp.InterpolateInplaceM(self.get_psf_model()[0], im, (fl['x'] - x0).astype('f4'), (fl['y'] - y0).astype('f4'), z_.astype('f4'), A.astype('f4'), roiSize, dx, dy, dz)
        else:
            for x_offset, z_offset, spec_chan, chan in zip(ChanXOffsets, ChanZOffsets, ChanSpecs, range(len(ChanSpecs))):
                z_ = np.clip(z - fl['z'] + z_offset, -maxz, maxz).astype('f')
                roiSize = np.minimum(8 + np.abs(z_) * (2.5 / dx), 140).astype('i')
                cInterp.InterpolateInplaceM(self.get_psf_model(chan)[0], im, (fl['x'] - x0 + x_offset).astype('f4'), (fl['y'] - y0).astype('f4'),
                                             z_.astype('f4'), (A * fl['spec'][:, spec_chan]).astype('f4'), roiSize, dx, dy, dz)

    def simPalmImFI(self, X, Y, z, fluors, intTime=.1, numSubSteps=10, laserPowers=[.1, 1], position=[0, 0, 0],
                     illuminationFunction='ConstIllum', ChanXOffsets=[0, ], ChanZOffsets=[0, ], ChanSpecs=None, im=None):
        if self.get_psf_model()[0] is None:
            self.genTheoreticalModel()

        if im is None:
            im = np.zeros((len(X), len(Y)), 'f')

        if fluors is None:
            return im

        # illuminationFunction may be a registry name (looked up in our own
        # illumination_functions dict) or an illumination-function instance
        if isinstance(illuminationFunction, str):
            illum_fcn = self.illumination_functions[illuminationFunction]
        else:
            illum_fcn = illuminationFunction

        A = np.zeros(len(fluors.fl), 'f')

        for n in range(numSubSteps):
            A += fluors.illuminate(laserPowers, intTime / numSubSteps, position=position, illuminationFunction=illum_fcn)

        dx = X[1] - X[0]
        dy = Y[1] - Y[0]
        dz = self.dz

        maxz = self.dz * (self.get_psf_model()[0].shape[2] / 2 - 1)

        x0 = X[0]
        y0 = Y[0]

        m = A > .1

        fl = fluors.fl[m]
        A2 = A[m]

        nCPUs = int(min(multiprocessing.cpu_count(), len(A2)))

        if nCPUs > 0:
            threads = [threading.Thread(target=self._render_fluor_subset,
                                         args=(im, fl[i::nCPUs], A2[i::nCPUs], x0, y0, z, dx, dy, dz, maxz, ChanXOffsets, ChanZOffsets, ChanSpecs))
                       for i in range(nCPUs)]

            for p in threads:
                p.start()

            for p in threads:
                p.join()

        return im
