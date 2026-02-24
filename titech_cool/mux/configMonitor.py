import json
import os
import datetime
from influxdb import InfluxDBClient

# insert config to influxDB
def upload(dir):
    scanContents = os.listdir(dir)
    # serch config file
    configfiles = [s for s in scanContents if ".json.befor" in s]
    for configfile in configfiles:
        # detect chip number
        if "Chip1" in configfile:
            chipNum = 1
        elif "Chip2" in configfile:
            chipNum = 2 
        elif "Chip3" in configfile:
            chipNum = 3 
        elif "Chip4" in configfile:
            chipNum = 4
        else:
            continue
        # read chipconfig json
        chipConf = json.load(open(dir+ "/"  + configfiles[0]))
        MonitorEnable = chipConf["RD53B"]["GlobalConfig"]["MonitorEnable"]
        MonitorV =  chipConf["RD53B"]["GlobalConfig"]["MonitorV"]
        # get timestamp
        timestamp = os.stat(dir).st_mtime
        date = datetime.datetime.fromtimestamp(timestamp)
        print(date)
        """ if you want to upload localdb, activate this part.
        try:
            client = InfluxDBClient( host = "atlastit01.kek.jp", port = 8086, database = "dcsDB" )
            data = [{
                "measurement":"TitechCoolingBox02",
                "fields":{
                    "MUX MonitorEnable ch"+str(chipNum) :MonitorEnable,"MUX MonitorV ch"+str(chipNum):MonitorV,
                }}]
            #res = client.write_points(data)
            print(data)
            client.close()
        except:
           print("cannot record date in influxDB")
        """

# read config
confFile = "/home/atlasj/titech_cool/mux/configMonitor.json"
config = json.load(open(confFile,"r"))

# read scan result directory list
scanDir = []
scanNum = []
for d in config["resulutDir"]: 
    for scan in os.listdir(d):
        if scan[0:6].isdigit():
            scanDir.append(d+"/"+scan)
            scanNum.append(int(scan[0:6]))
print(scanDir)
print(scanNum)

# detect unchecked result
for index,scan  in enumerate(scanNum):
    if scan > int(config["lastScanNumber"]):
        upload(scanDir[index])

# save
config["lastScanNumber"] = scan
with open(confFile,"w") as file:
    json.dump(config,file,indent=4)


