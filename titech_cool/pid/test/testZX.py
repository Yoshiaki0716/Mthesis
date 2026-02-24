#!/usr/bin/python3

from PeltierControl import TakasagoZXControl
import time

peltCtrl = TakasagoZXControl( { "port":"/dev/ttyUSB2", "address":1, "baud":9600 } )

peltCtrl.setOn()
peltCtrl.setOff()
print(peltCtrl.getID())
print(peltCtrl.getStatus())
peltCtrl.resetAlarms()
peltCtrl.setV(5.0)
peltCtrl.setOn()
time.sleep(3)
print(peltCtrl.getV())
print(peltCtrl.getI())
time.sleep(1)
peltCtrl.setOff()

