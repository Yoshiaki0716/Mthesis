#!/usr/bin/python3

from PeltierControl import TakasagoKXControl
import time

peltCtrl = TakasagoKXControl( { "port":"/dev/ttyUSB_TakasagoKX", "address":1, "baud":9600 } )

peltCtrl.setOn()
peltCtrl.setOff()
print(peltCtrl.getSetting())
peltCtrl.resetAlarms()
peltCtrl.setV(0.05)
print(peltCtrl.getSetting())
peltCtrl.setOn()
time.sleep(3)
print(peltCtrl.getV())
print(peltCtrl.getI())
print(peltCtrl.getStatus())
time.sleep(1)
peltCtrl.setV(0.0)
time.sleep(1)
print(peltCtrl.getSetting())
peltCtrl.setOff()
time.sleep(1)
print(peltCtrl.getSetting())

