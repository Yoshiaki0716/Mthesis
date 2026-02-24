#!/bin/python3

import json

with open('/tmp/watch.json') as f:
    j = json.load( f )
    values = j["results"][0]['series'][0]['values']
    columns = j["results"][0]['series'][0]['columns']
    for c,v in zip( columns, values[0] ):
        if c == 'Temperature (ch.target) [C]':
            setT = v
        if c == 'Set. Voltage (ch.Peltier) [V]':
            setV = v
    print( ' -t {} -V {}'.format(setT, setV) )
        
