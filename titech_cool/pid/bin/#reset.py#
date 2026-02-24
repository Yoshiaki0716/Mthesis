#!/usr/bin/env python3

import time
import sys
import optparse
import math

from PeltierControl import TexioPFRControl
from ArduinoSingleBoxDcsController import *
from ArduinoNTCReader import *


parser = optparse.OptionParser()
parser.add_option( '--zero', action = 'store_true', default = False, dest = 'resetPeltierVolt' )
options, remainder = parser.parse_args( sys.argv )

try:
    print('=========================\n\nArduinoSingleBoxDcsController Check\n\n');
    ctrl = ArduinoSingleBoxDcsController( { "port" : "/dev/ttyUSB_ArduinoSingleBoxController", "baud" : 9600 } )
    ctrl.reset()
    ctrl.readAll()
    ctrl.print()
    if math.isnan( ctrl.m_rh["carrier"] ):
        raise Exception('SHT85 RH readout is irregular!')
    print('==> OK')
except Exception as e:
    print('Exception raised in checking ArduinoDcsController: ', e)
    sys.exit(1)


for j in range(6):
    try:
        print('=========================\n\nTexio PFR\n\n');
        peltCtrl = TexioPFRControl( { "port":"/dev/ttyUSB_TexioPFR", "address":1, "baud":9600 } )
        peltCtrl.resetAlarms()
        time.sleep(1)
        
        if options.resetPeltierVolt:
            peltCtrl.setV( 0 )
            
        time.sleep(1)
        print( peltCtrl.getSetting() )
        print('==> OK')
        break
    except Exception as e:
        print('Exception raised in checking Texio PFR: ', e)
        if j == 2 :
            sys.exit(1)
        else:
            time.sleep(10)
            continue
            


print('\n\n==> All devices are green!')

sys.exit(0)
