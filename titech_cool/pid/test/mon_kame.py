#!/usr/bin/python3

from MCP3208 import *
from MAX31855 import *
from datetime import datetime

tMeters = {}

thermocouple = ThermocoupleReader( { 'pinout' : { 'ss' : 15 }, 'channel' : 0, 'calib': [ 0.0006892, -0.05584, 1.99 ] } )
tMeters.update( { 'peltier':thermocouple } )

channels = { "chiller" : 0, "chuck" : 2, "sink" : 5 }

for name,channel in channels.items():
    
    ntc = NtcReader( { "pinout"  : { "ss" : 31, "clk" : 33, "miso" : 35, "mosi" : 37 },
		       "channel" : channel,
                       "calib"   : { "R25" : 100, "Rext" : 470, "B" : 4281, "T0" : 298.81 } } )
    tMeters.update( { name : ntc } )




print( "# {} {} {} {} {} {}".format( "unixtime", "date",
                                     "peltier",
                                     "sink",
                                     "chuck",
                                     "chiller" ) )

while True:
    now = datetime.now()
    now_ts = now.timestamp()
    
    print( "{} {} {:.3f} {:.3f} {:.3f} {:.3f}".format( now_ts, now,
                                                       tMeters["peltier"].getT(),
                                                       tMeters["sink"]   .getT(),
                                                       tMeters["chuck"]  .getT(),
                                                       tMeters["chiller"].getT() ) )
    

    
