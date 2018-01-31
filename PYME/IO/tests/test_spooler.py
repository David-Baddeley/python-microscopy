import subprocess
from PYME.IO import clusterIO
from PYME.Acquire import HTTPSpooler
from PYME.IO import testClusterSpooling
import tempfile
import os
import shutil

from PYME.IO.clusterExport import ImageFrameSource, MDSource
from PYME.IO import MetaDataHandler
#import unittest
import time
import sys

procs = []
tmp_root = None

def setup_module():
    global proc, tmp_root
    tmp_root = os.path.join(tempfile.gettempdir(), 'PYMEDataServer_TEST')
    if not os.path.exists(tmp_root):
        os.makedirs(tmp_root)
    port_start = 8100
    for i in range(10):
        proc = subprocess.Popen('PYMEDataServer -r %s -f TEST -t -p %d --timeout-test=0' % (tmp_root, port_start + i), stderr= sys.stderr, shell=True)
        procs.append(proc)
        
    time.sleep(5)


def teardown_module():
    global proc, tmp_root
    #proc.send_signal(1)
    #time.sleep(1)
    for proc in procs:
        proc.kill()
    
    shutil.rmtree(tmp_root)
    
    
def test_spooler():
    ts = testClusterSpooling.TestSpooler(testFrameSize=[1024,256], serverfilter='TEST')
    ts.run(nFrames=2000)
    

from PYME.util import fProfile
    
if __name__ == '__main__':
    prof = fProfile.thread_profiler()
    setup_module()
    try:
        PROFILE = False
    
        if False:
            prof.profileOn('.*PYME.*|.*requests.*|.*socket.*|.*httplib.*', '/Users/david/spool_prof.txt')
            PROFILE = True
            
        test_spooler()

        if PROFILE:
            prof.profileOff()
        time.sleep(5)
    finally:
        teardown_module()