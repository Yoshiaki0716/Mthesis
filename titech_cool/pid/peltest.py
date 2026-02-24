from DeviceInvoker import *
from Logger import *
import json
import traceback
import time

with open("/home/admin/titech_cool/pid/config/RepicBox10.json") as f:
    config = json.load(f)
    devices = invokeDevices(config['devices'])

    peltCtrl = devices['peltierControl']
    peltCtrl.setOn()
    peltCtrl.setOff()
    peltCtrl.resetAlarms()

    peltCtrl.setOn()
    time.sleep(3)
    try:
        if abs( options.Vinit ) < 10.0 :
            if options.Vinit < 0.0:
                peltCtrl.setV( 0.0 )
                time.sleep(2)
                dcsCtrl.pelPolInverted()
            peltCtrl.setV( abs(options.Vinit) )
            pid.V     = options.Vinit
            pid.Vprev = options.Vinit
    except:
        print(traceback.format_exc())
       # logger.fatal( 'could not communicate with the peltier control power supply device. terminating' )
        print( 'could not communicate with the peltier control power supply device. terminating' )
        exit(1)

