#!/usr/bin/env python3

import time
import sys
import optparse

from PeltierControl import TexioPFRControl
from ArduinoSingleBoxDcsController import *
from ArduinoNTCReader import *


parser = optparse.OptionParser()
parser.add_option( '--zero', action = 'store_true', default = False, dest = 'resetPeltierVolt' )
options, remainder = parser.parse_args( sys.argv )

try:
    print('=========================\n\nArduinoSingleBoxDcsController Check\n\n');
    ctrl = ArduinoSingleBoxDcsController( { "port" : "/dev/ttyUSB_ArduinoSingleBoxController", "baud" : 9600 } )
    ctrl.readAll()
    ctrl.print()
    print('==> OK')
except Exception as e:
    print('Exception raised in checking ArduinoDcsController: ', e)
    sys.exit(1)

