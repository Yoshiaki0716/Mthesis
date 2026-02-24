import time
import math
from datetime import datetime
from InfluxAccess import *
import sys
import socket

if __name__ == '__main__':

    hostname = socket.gethostname()
    boxName = "CoolingBox{:02d}".format( int(hostname[-1]) )

    config_temp = { "host" : '192.168.100.104',
                    "port" : 8086,
                    'database': "REPIC",
                    "measurement": boxName,
                    "field" : "Temperature (ch.Module) [C]" }

    config_rh = { "host" : '192.168.100.104',
                  "port" : 8086,
                  'database': "REPIC",
                  "measurement": boxName,
                  "field" : "Humidity (ch.Carrier in) [%H]" }

    client_temp = InfluxAccess( config_temp )
    client_rh = InfluxAccess( config_rh )
    print( client_temp.getLast() )
    print( client_rh.getLast() )

