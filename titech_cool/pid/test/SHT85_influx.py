import time
import math
from datetime import datetime
from ArduinoSHT85Control import *

from influxdb import InfluxDBClient


client = InfluxDBClient( host = '192.168.10.128', port = 8086 )
client.switch_database( 'pidDB' )

sht85 = ArduinoSHT85Control( { "port" : "/dev/ttyACM2", "baud" : 9600 } )

now = datetime.now()
now_ts = now.timestamp()

temp = round(sht85.getT(), 3)
rh   = round(sht85.getRH(), 3)
dp   = round( 5.64*math.sqrt(rh) + 0.861*temp -  46.4, 3)
data = [ {'fields': { 'SHT85 Temp' :temp, 'SHT85 RH' : rh, 'SHT85 DP' : dp }, 'measurement': 'EnvMon' } ]

print( data )
res = client.write_points( data )
