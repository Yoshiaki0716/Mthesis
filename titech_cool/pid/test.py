from PidMonitor import *
import json
with open("/home/admin/titech_cool/pid/config/RepicBox10.json") as f:
    config = json.load(f)
    mon = PidMonitor( config['pidMonitor']['maxHistory'] )
