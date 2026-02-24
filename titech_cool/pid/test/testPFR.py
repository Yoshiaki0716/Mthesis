#!/usr/bin/env python3

from PeltierControl import TexioPFRControl
from ArduinoSingleBoxDcsController import ArduinoSingleBoxDcsController
import time
import matplotlib.pyplot as plt

peltCtrl = TexioPFRControl( { "port":"/dev/ttyUSB_TexioPFR", "address":1, "baud":9600 } )
ctrl = ArduinoSingleBoxDcsController( { "port" : "/dev/ttyUSB_ArduinoSingleBoxController", "baud" : 9600, 'useCalib':True } )

peltCtrl.idn()

def tryCommand( command, outType = str ):
    print( f'{command}: ', peltCtrl.command(command, outType) )

peltCtrl.setV(0.0)

print('------------------------------------')
tryCommand( 'output:protection:clear', None )
tryCommand( 'output:prot:trip?', bool )
tryCommand( 'volt:lim:auto OFF', None )
tryCommand( 'curr:lim:auto OFF', None )
tryCommand( 'output:mode?', int )
tryCommand( 'output?', bool )
tryCommand( 'volt?', float )
tryCommand( 'volt? max', float )
tryCommand( 'volt:lim:auto?', bool )
tryCommand( 'volt:prot?', float )
tryCommand( 'volt:prot? min', float )
tryCommand( 'curr?', float )
tryCommand( 'curr? max', float )
tryCommand( 'curr:lim:auto?', bool )
tryCommand( 'curr:prot?', bool )
tryCommand( 'curr:prot? min', float )
tryCommand( 'meas:current?', float )
tryCommand( 'meas:volt?', float )
print('------------------------------------')

IV_array = { "V" : [], "I" : [] }


peltCtrl.setOff()
print( 'setI = ', peltCtrl.setI(0.0) )
print( 'setV = ', peltCtrl.setV(0.0) )
print(peltCtrl.getSetting())
peltCtrl.resetAlarms()
time.sleep(1)
print('tripped = ', peltCtrl.isTripped())
print( 'setI = ', peltCtrl.setI(10.0) )
print( 'setV = ', peltCtrl.setV(1.0) )
time.sleep(1)

ctrl.relayOff( 4 );
ctrl.relayOff( 5 );

peltCtrl.setOn()
print( '==> set ON' )
for i in [0.5, 1.0, 2.0, 4.0, 6.5, 9.0]:
    print( 'setV = ', peltCtrl.setV(i) )
    time.sleep(1)
    print('status = ', peltCtrl.getStatus())
    print('tripped = ', peltCtrl.isTripped())
    print('V status = ', peltCtrl.getV(), ' [V]')
    print('I status = ', peltCtrl.getI(), ' [A]')
    
    IV_array["V"].append( peltCtrl.getV().get("meas") )
    IV_array["I"].append( peltCtrl.getI().get("meas") )
    print('------------')
    
peltCtrl.setOff()
print( '==> set OFF' )
time.sleep(2)

print( 'setV = ', peltCtrl.setV(0.0) )

ctrl.relayOn( 4 );
ctrl.relayOn( 5 );

peltCtrl.setOn()
print( '==> set ON' )
for i in [0.5, 1.0, 2.0, 4.0, 6.5, 9.0]:
    print( 'setV = ', peltCtrl.setV(i) )
    time.sleep(1)
    print('status = ', peltCtrl.getStatus())
    print('tripped = ', peltCtrl.isTripped())
    print('V status = ', peltCtrl.getV(), ' [V]')
    print('I status = ', peltCtrl.getI(), ' [A]')
    
    IV_array["V"].append( -1.0*peltCtrl.getV().get("meas") )
    IV_array["I"].append( -1.0*peltCtrl.getI().get("meas") )
    print('------------')
    
print( 'setV = ', peltCtrl.setV(0.0) )
peltCtrl.setOff()
print( '==> set OFF' )
time.sleep(2)

ctrl.relayOff( 4 );
ctrl.relayOff( 5 );

print( 'setV = ', peltCtrl.setV(0.0) )

print('status = ', peltCtrl.getStatus())
print('tripped = ', peltCtrl.isTripped())
print('V status = ', peltCtrl.getV(), ' [V]')
print('I status = ', peltCtrl.getI(), ' [A]')
print('setting: ', peltCtrl.getSetting())


plt.plot( IV_array["V"], IV_array["I"], marker='o', linestyle='None' )
plt.xlabel( "Peltier Voltage [V]" )
plt.ylabel( "Peltier Current [A]" )
plt.axis( ( -10, 10, -10, 10 ) )

plt.show()
