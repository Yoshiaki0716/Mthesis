#!/usr/bin/env python3

import time
from datetime import datetime
from ArduinoSHT85Control import *

sht85 = ArduinoSHT85Control( { "port" : "/dev/ttyACM0", "baud" : 9600 } )

while True:
    now = datetime.now()
    now_ts = now.timestamp()
    
    t = sht85.getT()
    rh = sht85.getRH()
    
    with open('/home/admin/test_REPIC_envlog.txt', 'a') as f:
        f.write( "{} {} {:.3f} {:.3f}\n".format( now_ts, now, t, rh ) )
        
    print( "{} {} {:.3f} {:.3f}\n".format( now_ts, now, t, rh ) )
    time.sleep(60)
