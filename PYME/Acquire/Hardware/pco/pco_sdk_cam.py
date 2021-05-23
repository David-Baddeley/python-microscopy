# -*- coding: utf-8 -*-

"""
Created on Sun May 16 2021

@author: zacsimile
"""

from PYME.Acquire.Hardware.Camera import Camera
from PYME.Acquire.Hardware.pco import pco_sdk
from PYME.Acquire import eventLog

import numpy as np
import logging
logger = logging.getLogger(__name__)

import ctypes
import ctypes.wintypes

import queue
import threading
import time

k32_dll = ctypes.windll.kernel32  # lets us use the recommended WaitForSingleObject call (see pco.sdk)
                                  # instead of the not-recommended-for-polling pco_sdk.get_buffer_status()

timebase = {pco_sdk.PCO_TIMEBASE_NS : 1e-9, 
            pco_sdk.PCO_TIMEBASE_US : 1e-6, 
            pco_sdk.PCO_TIMEBASE_MS : 1e-3}  # Conversions 

class camReg(object):
    """
    Keep track of the number of cameras initialised so we can initialise and
    finalise the library.
    """
    numCameras = -1

    @classmethod
    def regCamera(cls):
        if cls.numCameras == -1:
            pco_sdk.reset_lib()

        cls.numCameras += 1

    @classmethod
    def unregCamera(cls):
        cls.numCameras -= 1
        if cls.numCameras == 0:
            # There appears to be no pco.sdk uninitialization
            pass

camReg.regCamera()  # initialize/reset the sdk

MAX_BUFFERS = 16
MAX_TIMEOUTS = 10

class PcoSdkCam(Camera):
    def __init__(self, camNum, debuglevel='off'):
        Camera.__init__(self)
        self._initalized = False
        self.noiseProps = None

    def Init(self):
        self._handle = pco_sdk.open_camera()
        camReg.regCamera()
        self._desc = pco_sdk.get_camera_description(self._handle)

        if pco_sdk.get_recording_state(self._handle) == pco_sdk.PCO_CAMERA_RUNNING:
            pco_sdk.set_recording_state(self._handle, pco_sdk.PCO_CAMERA_STOPPED)
        pco_sdk.reset_settings_to_default(self._handle)
        _, err, _ = pco_sdk.get_camera_health_status(self._handle)
        if err != 0:
            self.Shutdown()
            raise pco_sdk.PcoSdkException(f"Camera shutdown with error status {err}.")

        self._bits_per_pixel = 16
        self._integ_time = 0
        self._delay_time = 0
        self._electr_temp = 0
        self._ccd_temp = 0
        self.__default_n_buffers = 16  # max of 16
        self._n_buffers = 0
        self._buf_event = None
        self._buf_addr = None
        self._recording = False
        self._roi = None
        self._timeout = 1000  # ms
        self._n_buffered = 0
        self._n_queued = 0
        self._binning_x = 0
        self._binning_y = 0
        self._n_timeouts = 0
        self._buffers_to_queue = queue.Queue()
        self._queued_buffers = queue.Queue()
        self._full_buffers = queue.Queue()
        self.SetROI(1, 1, self.GetCCDWidth(), self.GetCCDHeight())
        self.SetIntegTime(0.025)
        self.SetAcquisitionMode(self.MODE_CONTINUOUS)
        self._cam_type = pco_sdk.get_camera_type(self._handle)
        self.SetHotPixelCorrectionMode(pco_sdk.PCO_HOTPIXELCORRECTION_OFF)

        self._polling = True
        self._poll_thread = threading.Thread(target=self._poll_loop)
        self._poll_thread.start()

        self._initalized = True

    @property
    def noise_properties(self):
        return self.noiseProps

    def ExpReady(self):
        if not self._recording:
            return False

        return (self._n_buffered > 0)

    def _poll_loop(self):
        while self._polling:
            if self._recording:
                while not self._buffers_to_queue.empty():
                    lx, ly = self.GetPicWidth(), self.GetPicHeight()
                    i = self._buffers_to_queue.get()
                    pco_sdk.add_buffer_ex(self._handle, 0, i, lx, ly, self._bits_per_pixel)
                    self._queued_buffers.put(i)
                    self._n_queued += 1
                
                if self._n_queued > 0:
                    _curr_buf = self._queued_buffers.get()
                    self._n_queued -= 1
                    wait_status = k32_dll.WaitForSingleObject(self._buf_event[_curr_buf], self._timeout)
                    if wait_status:
                        self._n_timeouts += 1
                        if self._n_timeouts >= MAX_TIMEOUTS:
                            raise TimeoutError(f"Waited too long for buffer ({self._timeout} ms).")
                    k32_dll.ResetEvent(self._buf_event[_curr_buf])
                    self._full_buffers.put(_curr_buf)
                    self._n_buffered += 1
                else:
                    time.sleep(0.01)

                time.sleep(0.0001)
            else:
                time.sleep(0.01)

    def CamReady(self):
        return self._initalized

    def ExtractColor(self, chSlice, mode):
        if self._recording:
            _curr_buf = self._full_buffers.get()
            # update buffer index
            self._n_buffered -= 1

            # # Make sure this buffer works
            # _, status_drv = pco_sdk.get_buffer_status(self._handle, _curr_buf)
            # if status_drv:
            #     raise pco_sdk.PcoSdkException(f"Error {status_drv} during buffer check.")

            # copy image from _buf_addr
            ctypes.cdll.msvcrt.memcpy(chSlice.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16)),
                    self._buf_addr[_curr_buf], 
                    chSlice.nbytes)

            # recycle the buffer
            self._buffers_to_queue.put(_curr_buf)
        
    def GetName(self):
        return pco_sdk.get_camera_name(self._handle)

    def GetHeadModel(self):
        return pco_sdk.PCO_CAMERA_TYPES.get(self._cam_type.wCamType)

    def GetSerialNumber(self):
        return str(self._cam_type.dwSerialNumber)

    def SetIntegTime(self, time):
        lb = float(self._desc.dwMinExposDESC)*1e-9
        ub = float(self._desc.dwMaxExposDESC)*1e-3
        step = float(self._desc.dwMinExposStepDESC)*1e-9

        # Round to a multiple of time step
        time = np.round(time/step)*step

        # Don't let this go out of bounds
        time = np.clip(time, lb, ub)

        pco_sdk.set_delay_exposure_time(self._handle, int(self._delay_time*1e3), int(time*1e3), 
                                        pco_sdk.PCO_TIMEBASE_MS, pco_sdk.PCO_TIMEBASE_MS)
        _, exposure, _, timebase_exposure = pco_sdk.get_delay_exposure_time(self._handle)
        self._integ_time = exposure*timebase.get(timebase_exposure)

    def GetIntegTime(self):
        return self._integ_time

    def GetCycleTime(self):
        return (self._integ_time + self._delay_time)

    def GetCCDWidth(self):
        return self._desc.wMaxHorzResStdDESC

    def GetCCDHeight(self):
        return self._desc.wMaxVertResStdDESC

    def GetPicWidth(self):
        return self.GetROI()[2] - self.GetROI()[0] + 1

    def GetPicHeight(self):
        return self.GetROI()[3] - self.GetROI()[1] + 1

    def SetHorizontalBin(self, value):
        b_max, step_type = self._desc.wMaxBinHorzDESC, self._desc.wBinHorzSteppingDESC

        if step_type == 0:
            # binary binning, round to powers of two
            value = 1<<(value-1).bit_length()
        else:
            value = int(np.round(value)) # integer

        value = np.clip(value, 1, b_max)

        pco_sdk.set_binning(self._handle, value, self.GetVerticalBin())
        self._binning_x, _ = pco_sdk.get_binning(self._handle)

    def GetHorizontalBin(self):
        return self._binning_x

    def SetVerticalBin(self, value):
        b_max, step_type = self._desc.wMaxBinVertDESC, self._desc.wBinVertSteppingDESC

        if step_type == 0:
            # binary binning, round to powers of two
            value = 1<<(value-1).bit_length()
        else:
            value = int(np.round(value)) 

        value = np.clip(value, 1, b_max)

        pco_sdk.set_binning(self._handle, self.GetHorizontalBin(), value)
        _, self._binning_y = pco_sdk.get_binning(self._handle)

    def GetVerticalBin(self):
        return self._binning_y

    def GetSupportedBinning(self):
        import itertools

        bx_max, bx_step_type = self._desc.wMaxBinHorzDESC, self._desc.wBinHorzSteppingDESC
        by_max, by_step_type = self._desc.wMaxBinVertDESC, self._desc.wBinVertSteppingDESC

        if bx_step_type == 0:
            # binary step type
            x = [2**j for j in np.arange((bx_max).bit_length())]
        else:
            x = [j for j in np.arange(bx_max)]

        if by_step_type == 0:
            y = [2**j for j in np.arange((by_max).bit_length())]
        else:
            y = [j for j in np.arange(by_max)]

        return list(itertools.product(x,y))

    def SetROI(self, x0, y0, x1, y1):
        # Stepping (n pixels)
        dx = self._desc.wRoiHorStepsDESC
        dy = self._desc.wRoiVertStepsDESC

        # Chip size
        lx_max = self.GetCCDWidth()
        ly_max = self.GetCCDHeight()

        # Minimum ROI size
        lx_min = self._desc.wMinSizeHorzDESC
        ly_min = self._desc.wMinSizeVertDESC

        # Make sure bounds are ordered
        if x1 < x0:
            x0, x1 = x1, x0
        if y1 < y0:
            y0, y1 = y1, y0

        # Don't let the ROI go out of bounds
        # See pco.sdk manual chapter 3: IMAGE AREA SELECTION (ROI)
        x0 = np.clip(x0, 1, lx_max-dx+1)
        y0 = np.clip(y0, 1, ly_max-dy+1)
        x1 = np.clip(x1, 1+dx, lx_max)
        y1 = np.clip(y1, 1+dy, ly_max)

        # Don't let us choose too small an ROI
        if (x1-x0+1) < lx_min:
            logger.debug('Selected ROI width is too small, automatically adjusting to {}.'.format(lx_min))
            guess_pos = x0+lx_min
            # Deal with boundaries
            if guess_pos <= lx_max:
                x1 = guess_pos
            else:
                x0 = x1-lx_min-1
        if (y1-y0+1) < ly_min:
            logger.debug('Selected ROI height is too small, automatically adjusting to {}.'.format(ly_min))
            guess_pos = y0+ly_min
            if guess_pos <= ly_max:
                y1 = guess_pos
            else:
                y0 = y1-ly_min-1

        # Round to a multiple of dx, dy
        # TODO: Why do I need the 1 correction only on x0, y0???
        x0 = 1+int(np.floor((x0-1)/dx)*dx)
        y0 = 1+int(np.floor((y0-1)/dy)*dy)
        x1 = int(np.floor(x1/dx)*dx)
        y1 = int(np.floor(y1/dy)*dy)

        pco_sdk.set_roi(self._handle, x0, y0, x1, y1)
        self._roi = [x0, y0, x1, y1]

        logger.debug('ROI set: x0 %3.1f, y0 %3.1f, w %3.1f, h %3.1f' % (x0, y0, x1-x0+1, y1-y0+1))

    def GetROI(self):
        return self._roi
    
    def GetElectrTemp(self):
        return self._electr_temp

    def GetCCDTemp(self):
        return self._ccd_temp

    def GetCCDTempSetPoint(self):
        return self._ccd_temp_set_point

    def SetCCDTemp(self, temp):
        lb = self._desc.sMinCoolSetDESC
        ub = self._desc.sMaxCoolSetDESC
        temp = np.clip(temp, lb, ub)

        pco_sdk.set_cooling_setpoint_temperature(self._handle, temp)
        self._ccd_temp_set_point = temp
        self._get_temps()

    def _get_temps(self):
        # NOTE: temperature only gets probed when acquisition starts/stops (which can be fairly
        # irregularly - to the point of not being useful).
        # TODO - find a way to safely call this while the camera is running
        # FIXME: should _electr_tmp be power temperature rather than cam temp?
        self._ccd_temp, self._electr_temp, _ = pco_sdk.get_temperature(self._handle)

    def GetAcquisitionMode(self):
        return self._mode

    def SetAcquisitionMode(self, mode):
        if mode in [self.MODE_CONTINUOUS, self.MODE_SINGLE_SHOT]:
            self._mode = mode
        else:
            raise RuntimeError(f"Mode {mode} not supported")
        if mode == self.MODE_CONTINUOUS:
            self._n_buffers = self.__default_n_buffers
        elif mode == self.MODE_SINGLE_SHOT:
            self._n_buffers = 2

    def StartExposure(self):
        self.StopAq()
        
        lx, ly = self.GetPicWidth(), self.GetPicHeight()
        self._timeout = int(max(2*100*self.GetCycleTime(), 100))
        bufsize = lx*ly*ctypes.sizeof(ctypes.wintypes.WORD)  # how many words is this image worth?
        for i in np.arange(self._n_buffers):
            self._buf_event.append(ctypes.c_void_p(0))
            self._buf_addr.append(ctypes.POINTER(ctypes.wintypes.WORD)())
            pco_sdk.allocate_buffer(self._handle, -1, bufsize, 
                                    self._buf_addr[i], self._buf_event[i])
            self._buffers_to_queue.put(i)

        pco_sdk.set_image_parameters(self._handle, lx, ly, pco_sdk.PCO_IMAGEPARAMETERS_READ_WHILE_RECORDING)
        pco_sdk.arm_camera(self._handle)
        pco_sdk.set_recording_state(self._handle, pco_sdk.PCO_CAMERA_RUNNING)
        eventLog.logEvent('StartAq', '')
        self._recording = True

        return 0

    def StopAq(self):
        if self._recording:
            self._recording = False
            pco_sdk.set_recording_state(self._handle, pco_sdk.PCO_CAMERA_STOPPED)
            pco_sdk.cancel_images(self._handle)
            for i in np.arange(self._n_buffers):
                pco_sdk.free_buffer(self._handle, i)
            while not self._buffers_to_queue.empty():
                self._buffers_to_queue.get()
            while not self._queued_buffers.empty():
                self._queued_buffers.get()
            while not self._full_buffers.empty():
                self._full_buffers.get()
        self._n_buffered = 0
        self._n_queued = 0
        self._n_timeouts = 0
        self._buf_event = []
        self._buf_addr = []
        self._get_temps()

    def GetNumImsBuffered(self):
        return self._n_buffered

    def GetBufferSize(self):
        return self._n_buffers

    def SetBufferSize(self, n_buffers):
        if n_buffers > MAX_BUFFERS:
            logger.debug(f"{n_buffers} is greater than the maximum number of buffers, {MAX_BUFFERS}. Defaulting to {MAX_BUFFERS}.")
            n_buffers = MAX_BUFFERS
        self.__default_n_buffers = n_buffers
        self._n_buffers = n_buffers

    def GetFPS(self):
        if self.GetCycleTime() == 0:
            return 0
        return 1.0/self.GetCycleTime()

    def SetHotPixelCorrectionMode(self, mode):
        pco_sdk.set_hot_pixel_correction_mode(self._handle, mode)

    def SetBitsPerPixel(self, bits_per_pixel):
        self._bits_per_pixel = bits_per_pixel

    def Shutdown(self):
        self._polling = False
        pco_sdk.close_camera(self._handle)
        camReg.unregCamera()
