#!/usr/bin/python3

from MCP3208 import *

import sys

channel = int( sys.argv[1] )

ntc = NtcReader( { 
                   'pinout' : { 'ss': 31, 'clk': 37, 'miso': 35, 'mosi': 33, 'w': 1.e-3 },
                   'channel' : channel,
    'calib' : { 'R25' : 100, 'Rext' : 470, 'B':4281, 'T0':273.15+27.16 } } )

for i in range(20):
    print( 'channel = {}, T = {}'.format( channel, ntc.getT() ) )



    
