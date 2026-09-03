#!/usr/bin/python

##################
# cameranoise.py
#
# Copyright David Baddeley, 2009
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

import numpy as np

from PYME.misc import EMCCDTheory


def generate_camera_maps(size_x = 1024, size_y = 1024, seed=100, read_median=1.38, offset=100):
    """
    Generate camera maps for sCMOS simulation, using a constant random seed so that the maps are reproducible
    
    The use (and parameterization) of pareto distributions is designed to match the distribution of values observed in
    actual camera maps. Note that the pareto gives a somewhat better match than lognormal.
    
    Parameters
    ----------
    size_x
    size_y
    seed
    read_median
    offset

    Returns
    -------

    """
    
    np.random.seed(seed)
    
    #Variance
    s = 2.0
    #var = (np.random.lognormal(np.log(read_median), s, [size_x, size_y]))**2
    var = (read_median/(2**(1./s)) * (1 + np.random.pareto(s, [size_x, size_y]))) ** 2
    
    #the dark map has 3 components - a pareto distributed base distribution, a small ammount of Gaussian spread, and Gaussian distributed fixed pattern
    # line noise
    dark = offset + np.random.pareto(2.7, [size_x, size_y]) + np.random.normal(0, 1.8, [size_x, size_y]) + np.random.normal(0, 0.35, [size_x,])[:,None]
    
    flatfield = np.ones_like(dark)
    
    np.random.seed()
    return {'variance': var, 'dark':dark, 'flat' : flatfield}

class NoiseMaker:
    def __init__(self, QE=.8, electronsPerCount=27.32, readoutNoise=109.8, EMGain=0, background=0., floor=967, shutterOpen = True,
                 numGainElements=536, vbreakdown=6.6, temperature = -70., fast_read_approx=True):
        self.QE = QE
        self.ElectronsPerCount = electronsPerCount
        self.ReadoutNoise=readoutNoise
        self.EMGain=EMGain
        self.background = background
        self.ADOffset = floor
        self.NGainElements = numGainElements
        self.vbreakdown = vbreakdown
        self.temperature = temperature
        self.shutterOpen = shutterOpen
        
        self.approximate_read_noise = fast_read_approx #approximate readout noise
        
        self._ar_key = None
        self._ar_cache = None
        
    def _read_approx(self, im_shape):
        """
        Really dirty fast approximation to readout noise by indexing into a random location within a pre-calculated noise
        matrix. Note that this may result in undesired correlations in the read noise.
        
        Parameters
        ----------
        im_shape

        Returns
        -------

        """
        nEntries = int(np.prod(im_shape))
        ar_key = (nEntries, self.ADOffset, self.ReadoutNoise, self.ElectronsPerCount)
        
        if not self._ar_key == ar_key or self._ar_cache is None:
            self._ar_cache = self.ADOffset + (self.ReadoutNoise / self.ElectronsPerCount)*np.random.normal(size=2*nEntries)
            self._ar_key = ar_key
            
        offset = np.random.randint(0, nEntries)
        return self._ar_cache[offset:(offset+nEntries)].reshape(im_shape)

    def noisify(self, im):
        """Add noise to image using an EMCCD noise model
        
        Inputs
        ------
        
        im : NxM array of intensities (in photons)
        
        Outputs
        -------
        
        out: NxM array of simulated camera pixel intensities (in ADUs)
        
        """

        M = EMCCDTheory.M((80. + self.EMGain)/(255 + 80.), self.vbreakdown, self.temperature, self.NGainElements, 2.2)
        F2 = 1.0/EMCCDTheory.FSquared(M, self.NGainElements)

        if self.approximate_read_noise:
            o = self._read_approx(im.shape)
        else:
            o = self.ADOffset + (self.ReadoutNoise / self.ElectronsPerCount) * np.random.standard_normal(im.shape)
        
        if self.shutterOpen:
            o = o +  (M/(self.ElectronsPerCount*F2))*np.random.poisson((self.QE*F2)*(im + self.background))

        return o
        
    def getbg(self):
        M = EMCCDTheory.M((80. + self.EMGain)/(255 + 80.), self.vbreakdown, self.temperature, self.NGainElements, 2.2)
        F2 = 1.0/EMCCDTheory.FSquared(M, self.NGainElements)

        return self.ADOffset + M*(int(self.shutterOpen)*(0 + self.background)*self.QE*F2)/(self.ElectronsPerCount*F2) 
