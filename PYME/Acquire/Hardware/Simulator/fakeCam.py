#!/usr/bin/python

##################
# fakeCam.py
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

from PYME.simulation import image_gen
from PYME.simulation.cameranoise import generate_camera_maps, NoiseMaker
import scipy

from PYME.IO import MetaDataHandler
from PYME.Acquire import eventLog

import numpy as np

import threading
#import processing
import time

import ctypes
import sys

if sys.platform == 'win32':
    memcpy = ctypes.cdll.msvcrt.memcpy
elif sys.platform == 'darwin':
    memcpy = ctypes.CDLL('libSystem.dylib').memcpy
else: #linux
    memcpy = ctypes.CDLL('libc.so.6').memcpy

from PYME.Acquire.Hardware import ccdCalibrator

import logging
logger = logging.getLogger(__name__)

WELL_DEPTH= (2 << 15) -1


#calculate image in a separate thread to maintain GUI reponsiveness
class compThread(threading.Thread):
    def __init__(self,XVals, YVals,zPiezo, zOffset, fluors, noisemaker, laserPowers, intTime, contMode = True,
                 bufferlength=500, biplane = False, biplane_z = 500, xpiezo=None, ypiezo=None, illumFcn = 'ConstIllum', objects=None, image_generator=None):
        #TODO - Do we need to change the default buffer length. This shouldn't really be an issue as we pause the simulation the buffer starts to fill up.
        
        threading.Thread.__init__(self)
        self.XVals = XVals
        self.YVals = YVals
        self.fluors = fluors # type: PYME.Acquire.Hardware.Simulator.fluor.Fluorophores
        self.objects = objects #list of Fluorophores instances
        #self.zPos = zPos
        self.laserPowers = laserPowers
        self.intTime = intTime
        self.noiseMaker = noisemaker
        self.contMode = contMode
        self.bufferlength = bufferlength
        self.buffer = np.zeros((bufferlength, len(XVals), len(YVals)), 'uint16')
        self.bufferWritePos = 0
        self.bufferReadPos = 0
        self.numBufferedImages = 0

        assert(image_generator is not None)
        self.image_generator = image_generator

        self.biplane = biplane
        self.deltaZ = biplane_z

        self.zPiezo = zPiezo
        self.zOffset = zOffset

        self.xPiezo = xpiezo
        self.yPiezo = ypiezo
        self.illumFcn = illumFcn

        self.kill = False
        self.aqRunning = False
        self.stopAq = False
        self.startAq = False

    def setSplitterInfo(self, chan_z_offsets, chan_specs, chan_x_offsets=None):
        self._chan_z_offsets = chan_z_offsets
        self._chan_specs = chan_specs

        if chan_x_offsets:
            self._chan_x_offsets = chan_x_offsets
        else:
            nChans = len(chan_z_offsets)
            x_pixels = len(self.XVals)
            x_chan_pixels = x_pixels/nChans
            x_chan_size = (self.XVals[1] - self.XVals[0])*x_chan_pixels

            self._chan_x_offsets = [i*x_chan_size for i in range(nChans)]

    @property
    def ChanXOffsets(self):
        try:
            return getattr(self, '_chan_x_offsets')
        except AttributeError:
            if not self.fluors:
                return [0,]
            elif not self.biplane and not 'spec' in self.fluors.fl.dtype.fields.keys():
                return [0,]
            else:
                return [0, self.XVals[self.XVals.shape[0] / 2] - self.XVals[0]]

    @property
    def ChanZOffsets(self):
        try:
            return getattr(self, '_chan_z_offsets')
        except AttributeError:
            if not self.fluors:
                return [0, ]
            elif not self.biplane and not 'spec' in self.fluors.fl.dtype.fields.keys():
                return [0, ]
            else:
                return [0, self.deltaZ]

    @property
    def ChanSpecs(self):
        try:
            return getattr(self, '_chan_specs')
        except AttributeError:
            if not self.fluors:
                return None
            elif not 'spec' in self.fluors.fl.dtype.fields.keys():
                return None
            else:
                return [0,1]


    def run(self):
        while not self.kill:
            #self.frameLock.acquire()
            while ((not self.aqRunning) or (self.numBufferedImages > self.bufferlength/2.)) and (not self.kill) :
                time.sleep(.01)

            zPos = (self.zPiezo.effective_pos - self.zOffset)*1e3

            xp = 0
            yp = 0
            if not self.xPiezo is None:
                xp = (self.xPiezo.effective_pos - self.xPiezo.max_travel/2)*1e3

            if not self.xPiezo is None:
                yp = (self.yPiezo.effective_pos - self.yPiezo.max_travel/2)*1e3

            roi_bbox = (xp,yp, xp+self.XVals[-1], yp + self.YVals[-1])

            #print self.ChanSpecs, self.ChanXOffsets

            if self.objects is not None:
                r_i = np.zeros((len(self.XVals), len(self.YVals)), 'f')
                for obj in self.objects:
                    if obj.hit_test(roi_bbox):
                        self.image_generator.simulate_image(self.XVals, self.YVals,obj,
                                                                  laserPowers=self.laserPowers, intTime=self.intTime,
                                                                  position=[xp,yp,zPos], illuminationFunction=self.illumFcn,
                                                                  ChanXOffsets=self.ChanXOffsets, ChanZOffsets=self.ChanZOffsets,
                                                                  ChanSpecs=self.ChanSpecs, im=r_i)
            else:
                r_i = self.image_generator.simulate_image(self.XVals, self.YVals, self.fluors,
                                                                  laserPowers=self.laserPowers, intTime=self.intTime,
                                                                  position=[xp,yp,zPos], illuminationFunction=self.illumFcn,
                                                                  ChanXOffsets=self.ChanXOffsets, ChanZOffsets=self.ChanZOffsets,
                                                                  ChanSpecs=self.ChanSpecs)
                        
            r_i = r_i[:,:]
            _im = self.noiseMaker.noisify(r_i)
            self.im = np.clip(_im, 0, WELL_DEPTH).astype('uint16')

            self.buffer[self.bufferWritePos,:,:] = self.im
            self.bufferWritePos +=1
            if self.bufferWritePos >= self.bufferlength: #wrap around
                self.bufferWritePos = 0

            self.numBufferedImages = min(self.numBufferedImages +1, self.bufferlength)


            if not self.contMode:
                self.aqRunning = False

            if self.stopAq:
                self.aqRunning = False
                self.bufferWritePos = 0
                self.bufferReadPos = 0
                self.numBufferedImages = 0
                self.stopAq = False

            if self.startAq:
                self.aqRunning = True
                self.startAq = False

            #self.frameLock.release()

    def numFramesBuffered(self):
        return self.numBufferedImages

    def StartExp(self):
        self.bufferWritePos = 0
        self.bufferReadPos = 0
        self.numBufferedImages = 0
        self.aqRunning = True
        self.startAq = True
        #self.frameLock.release()

    def getIm(self):
        im = np.copy(self.buffer[self.bufferReadPos,:,:], order='F')
        self.numBufferedImages -= 1
        self.bufferReadPos +=1
        if self.bufferReadPos >= self.bufferlength: #wrap around
            self.bufferReadPos = 0

        return im

    def StopAq(self):
        self.stopAq = True
#        self.aqRunning = False
#        self.bufferWritePos = 0
#        self.bufferReadPos = 0
#        self.numBufferedImages = 0

        




from PYME.Acquire.Hardware.Camera import Camera
class FakeCamera(Camera):
    numpy_frames=1
    order= 'C'
    #MODE_CONTINUOUS=True
    #MODE_SINGLE_SHOT=False
    
    def __init__(self, XVals, YVals, noiseMaker, zPiezo, zOffset=50.0, fluors=None, laserPowers=[0,50], xpiezo=None, ypiezo=None, illumFcn = 'ConstIllum', pixel_size_nm=70.):
        if np.isscalar(XVals):
            self.SetSensorDimensions(XVals, YVals, pixel_size_nm, restart=False)
            self.pixel_size_nm = pixel_size_nm
        else:
            self.XVals = XVals
            self.YVals = YVals
    
            self.ROIx = (0,len(XVals))
            self.ROIy = (0,len(YVals))
            
            self.pixel_size_nm = XVals[1] - XVals[0]


        self.image_generator = image_gen.ImageGenerator()
        self.image_generator.set_pixelsize_nm(self.pixel_size_nm)

        self.zPiezo=zPiezo
        self.xPiezo = xpiezo
        self.yPiezo = ypiezo
        self.fluors=fluors
        self._objects=None
        self.noiseMaker=noiseMaker

        self._saturation_threshold = (2**16) - 1
        self.DefaultEMGain = 150
        self.preampGain = 1
        
        self.laserPowers=laserPowers
        self.illumFcn = illumFcn

        self.intTime=0.1
        self.zOffset = zOffset

        self.compT = None #thread which is currently being computed
        self._restart_compT()

        self._acquisition_mode = self.MODE_CONTINUOUS
        #self.contMode = True
        self.shutterOpen = True

        #let us work with andor dialog
        self.HorizShiftSpeeds = [[[10]]]
        self.vertShiftSpeeds = [1]
        self.fastestRecVSInd = 0
        self.frameTransferMode = False
        self.HSSpeed = 0
        self.VSSpeed = 0

        self.active = True

        #register as a provider of metadata
        MetaDataHandler.provideStartMetadata.append(self.GenStartMetadata)

        Camera.__init__(self)

    def setSplitterInfo(self, chan_z_offsets, chan_specs):
        self._chan_z_offsets = chan_z_offsets
        self._chan_specs = chan_specs


    def setFluors(self, fluors):
        self.fluors = fluors

        self._restart_compT()

    def set_objects(self, objs):
        self._objects = objs

        self._restart_compT()
        
    def SetSensorDimensions(self, x_size=256, y_size=256, pixel_size_nm=70., restart=True):
        self.XVals = pixel_size_nm*np.arange(0.0, float(x_size))
        self.YVals = pixel_size_nm * np.arange(0.0, float(y_size))
            
        self.ROIx = (0, len(self.XVals))
        self.ROIy = (0, len(self.YVals))
        
        if restart:
            self._restart_compT()

    def _preamp_mode_repr(self):
        return 'Preamp mode %d' % self.preampGain
    
    def GetSerialNumber(self):
        return 'FAKE-000'
    
    def SetIntegTime(self, iTime): 
        self.intTime=iTime#*1e-3
        self.compT.intTime = iTime#*1e-3
    def GetIntegTime(self): 
        return self.intTime
    
    def GetCCDWidth(self): 
        return len(self.XVals)
    def GetCCDHeight(self): 
        return len(self.YVals)
    
    def GetCCDTemp(self):
        return self.noiseMaker.temperature
    
    def GetPicWidth(self): 
        return self.ROIx[1] - self.ROIx[0]
    def GetPicHeight(self): 
        return self.ROIy[1] - self.ROIy[0]

    def SetROI(self, x1, y1, x2, y2):
        self.ROIx = (x1, x2)
        self.ROIy = (y1, y2)
        
        self._restart_compT()
        
    def _restart_compT(self):
        try:
            running = self.compT.aqRunning
            self.compT.kill = True
            while self.compT.is_alive():
                time.sleep(0.01)
                
        except AttributeError:
            running = False


        if self._objects is not None:
            self.compT = compThread(self.XVals[self.ROIx[0]:self.ROIx[1]], self.YVals[self.ROIy[0]:self.ROIy[1]],
                                self.zPiezo, self.zOffset, None, self.noiseMaker, laserPowers=self.laserPowers,
                                intTime=self.intTime, xpiezo=self.xPiezo, ypiezo=self.yPiezo, illumFcn=self.illumFcn, 
                                objects=self._objects, image_generator=self.image_generator)
        else:
            self.compT = compThread(self.XVals[self.ROIx[0]:self.ROIx[1]], self.YVals[self.ROIy[0]:self.ROIy[1]],
                                self.zPiezo, self.zOffset, self.fluors, self.noiseMaker, laserPowers=self.laserPowers,
                                intTime=self.intTime, xpiezo=self.xPiezo, ypiezo=self.yPiezo, illumFcn=self.illumFcn,
                                image_generator=self.image_generator)

        try:
            #vx = self.XVals[1] - self.XVals[0]
            chan_x_offsets = getattr(self, '_chan_x_offsets', None)
            print('chan_x_offsets:', chan_x_offsets)
            #except AttributeError:
            #    chan_x_offsets=None

            self.compT.setSplitterInfo(self._chan_z_offsets, self._chan_specs, chan_x_offsets=chan_x_offsets)
        except AttributeError:
            pass
        
        self.compT.start()
        if running:
            self.compT.StartExp()

        #self.compT.aqRunning = running
        
    def GetROI(self):
        return self.ROIx[0], self.ROIy[0], self.ROIx[1], self.ROIy[1]


    def Shutdown(self):
        self.compT.kill = True
        #pass

    def StartAq(self):
        self.compT.StartExp()
        #pass

    def StopAq(self):
        self.compT.StopAq()
        #pass

    def StartExposure(self):
        self._log_exposure_start()
        self.compT.StartExp()
        return 0


    def ExpReady(self):
        #return not self.compTCur.isAlive() #thread has finished -> a picture is available
        return self.compT.numFramesBuffered() > 0
 
    def ExtractColor(self, chSlice, mode): 
        try:
            d = self.compT.getIm()
            #print d.nbytes, chSlice.nbytes
            memcpy(chSlice.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
                   d.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)), chSlice.nbytes)
            #chSlice[:,:] = d #grab image from completed computation thread
            #self.compTOld = None #set computation thread to None such that we get an error if we try and obtain the same result twice
        except AttributeError:  # triggered if called with None
            logger.error("Grabbing problem: probably called with 'None' thread")
        
        
    def GetNumImsBuffered(self):
        return self.compT.numFramesBuffered()
    
    def GetBufferSize(self):
        return self.compT.bufferlength

    def GenStartMetadata(self, mdh):
        self.GetStatus()

        mdh.setEntry('Camera.Name', 'Simulated EM CCD Camera')

        mdh.setEntry('Camera.IntegrationTime', self.GetIntegTime())
        mdh.setEntry('Camera.CycleTime', self.GetIntegTime())
        mdh.setEntry('Camera.EMGain', self.GetEMGain())

        #mdh.setEntry('Camera.ROIPosX', self.GetROIX1())
        #mdh.setEntry('Camera.ROIPosY',  self.GetROIY1())
        
        x1, y1, x2, y2 = self.GetROI()
        mdh.setEntry('Camera.ROIOriginX', x1)
        mdh.setEntry('Camera.ROIOriginY', y1)
        mdh.setEntry('Camera.ROIWidth', x2 - x1)
        mdh.setEntry('Camera.ROIHeight', y2 - y1)
        #mdh.setEntry('Camera.StartCCDTemp',  self.GetCCDTemp())

        mdh.setEntry('Camera.ReadNoise', self.noiseMaker.ReadoutNoise)
        mdh.setEntry('Camera.NoiseFactor', 1.41)
        mdh.setEntry('Camera.ElectronsPerCount', self.noiseMaker.ElectronsPerCount)
        mdh.setEntry('Camera.ADOffset', np.mean(self.noiseMaker.ADOffset))

        #mdh.setEntry('Simulation.Fluorophores', self.fluors.fl)
        #mdh.setEntry('Simulation.LaserPowers', self.laserPowers)

        realEMGain = ccdCalibrator.getCalibratedCCDGain(self.GetEMGain(), self.GetCCDTempSetPoint())
        if not realEMGain is None:
            mdh.setEntry('Camera.TrueEMGain', realEMGain)
            
        if self.fluors and 'spec' in self.fluors.fl.dtype.fields.keys(): #set the splitter parameters
            mdh['Splitter.Channel0ROI'] = [0,0,128, 256]
            mdh['Splitter.Channel1ROI'] = [128,0,128, 256]
            mdh['Splitter.Flip'] = False

        chan_specs = getattr(self, '_chan_specs', None)
        if not chan_specs is None:
            nChans  = len(chan_specs)
            x_pixels = len(self.XVals)
            x_chan_pixels = x_pixels / nChans
            y_pixels = len(self.YVals)
            mdh['Multiview.NumROIs'] = nChans
            mdh['Multiview.ROISize'] =  [x_chan_pixels, y_pixels]
            mdh['Multiview.ChannelColor'] =  list(chan_specs)
            mdh['Splitter.Flip'] = False

            # write shift information (zero shifts)
            mdh['chroma.dx'] = '{"PYME.Analysis.points.twoColour.lin2Model": {"mx": 0, "my": 0, "x0": 0}}'
            mdh['chroma.dy'] = '{"PYME.Analysis.points.twoColour.lin2Model": {"mx": 0, "my": 0, "x0": 0}}'

            for i in range(nChans):
                mdh['Multiview.ROI%dOrigin' % i] = [i*x_chan_pixels, 0]
                mdh['Splitter.Channel%dROI' % i] = [i*x_chan_pixels, 0, x_chan_pixels, y_pixels]

        
            mdh['Simulator.ChanZOffsets'] = self._chan_z_offsets
            mdh['Simulator.ChanSpecs'] = self._chan_specs

    #functions to make us look more like andor camera
    def GetEMGain(self):
        return self.noiseMaker.EMGain

    def GetCCDTempSetPoint(self):
        return self.GetCCDTemp()

    def SetCCDTemp(self, temp):
        self.noiseMaker.temperature = temp
        #pass

    def SetEMGain(self, gain):
        self.noiseMaker.EMGain = gain
        #pass

    def GetAcquisitionMode(self):
        return self._acquisition_mode
    
    def SetAcquisitionMode(self, mode):
        self._acquisition_mode = mode
        self.compT.contMode = (mode == self.MODE_CONTINUOUS)

    def SetShutter(self, mode):
        self.shutterOpen = mode
        self.noiseMaker.shutterOpen = mode

    def GetBaselineClamp(self):
        return True


    def SetIlluminationFcn(self, illumFcn):
        self.illumFcn = illumFcn
        self.compT.illumFcn = illumFcn

    def __getattr__(self, name):
        if name in dir(self.noiseMaker):
            return self.noiseMaker.__dict__[name]
        else:  raise AttributeError(name)  # <<< DON'T FORGET THIS LINE !!
        
    def __del__(self):
        self.Shutdown()
        #self.compT.kill = True
