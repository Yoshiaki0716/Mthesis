#!/usr/bin/env python3

import time
import sys
import optparse

from InfluxAccess import *
import os
import json
import pprint
import datetime

pid_path=os.environ["PIDPATH"]
with open(f"{pid_path}/config/default.json") as f:
    default_config = json.load(f)

cfg_file = default_config.get("config")
with open(f"{pid_path}/{cfg_file}") as f:
    config = json.load(f)


fields = [ "isPidActivated (ch.sum)" ]

meas = config.get("id")


config.get("id")
ctrl = InfluxAccess( { "host" : '192.168.100.104', "port" : 8086, "database" : "REPIC",
                       "measurement" : meas, "field" : fields[0] } )

last=datetime.datetime.fromisoformat( ctrl.getLastTime() ).timestamp()
now=datetime.datetime.now().timestamp()

dt = now - last

print( dt )
