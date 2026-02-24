#!/usr/bin/env python3

from ArduinoSingleBoxDcsController import ArduinoSingleBoxDcsController
import os, sys

os.system("kill -9 $( ps ux | grep -i testread | grep -v grep | awk '{ print $2; }' )" )

cmd =  ' '.join( sys.argv[1:] )

ctrl = ArduinoSingleBoxDcsController( { "port" : "/dev/ttyUSB_ArduinoSingleBoxController", "baud" : 9600, 'useCalib':True } )

from PeltierControl import TexioPFRControl
try:
    peltCtrl = TexioPFRControl( { "port":"/dev/ttyUSB_TexioPFR", "address":1, "baud":9600 } )
except:
    peltCtrl = None


if sys.argv[1] == 'setV':
    V = float( sys.argv[2] )
    peltCtrl.setV( 0.0 )
    
    if V > 0:
        ctrl.pelPolNormal()
        peltCtrl.setV( abs(V) )
    else:
        ctrl.pelPolInverted()
        peltCtrl.setV( abs(V) )

    peltCtrl.setOn()
    
chmap = { 'heater':1, 'peltier':2, 'chiller':3, 'pelP' : 4, 'pelN' : 5, 'LV' : 6, 'HV' : 7, 'Lock' : 8 }

if sys.argv[1] == 'relayOff':
    ctrl.relayOff( chmap[sys.argv[2]] )

if sys.argv[1] == 'relayOn':
    ctrl.relayOn( chmap[sys.argv[2]] )

if sys.argv[1] == 'reset':
    ctrl.reset()

if sys.argv[1] == 'unlock':
    ctrl.unlock()

del ctrl
del peltCtrl

os.system("cd ${PIDPATH}; ./test/testRead.py >/dev/null &" )
