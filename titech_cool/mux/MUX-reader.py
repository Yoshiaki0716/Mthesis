import subprocess
import os
import time
import datetime
import json
import pandas as pd
from influxdb import DataFrameClient, InfluxDBClient

########### config
confFile = "/home/atlasj/titech_cool/mux/MUX-reader.json"
config = json.load(open(confFile,"r"))

YARR_DIR =config["YARR_DIR"] 
controller = config["controller"]
connectivity = config["connectivity"]
MUX_list=[int(i) for i in config["muxList"]]

###########


result={
    "Chip1":{},
    "Chip2":{},
    "Chip3":{},
    "Chip4":{}
}

os.chdir(YARR_DIR)
DFclient = DataFrameClient( host = "atlastit01.kek.jp", port = 8086, database = "dcsDB" )
client = InfluxDBClient( host = "atlastit01.kek.jp", port = 8086, database = "dcsDB") 
def upload_conf(MonitorEnable,MonitorV):
    data = [{
            "measurement":"TitechCoolingBox02",
            "fields":{
                "MUX MonitorEnable ch1":MonitorEnable,"MUX MonitorV ch1":MonitorV,
                "MUX MonitorEnable ch2":MonitorEnable,"MUX MonitorV ch2":MonitorV,
                "MUX MonitorEnable ch3":MonitorEnable,"MUX MonitorV ch3":MonitorV,
                "MUX MonitorEnable ch4":MonitorEnable,"MUX MonitorV ch4":MonitorV}
            }]
    res = client.write_points(data)
  
  
  
def getMUX(configChangeDate):
    flag = True
    while flag:
        MUX_1 = DFclient.query("SELECT * FROM TitechCoolingBox02 WHERE 'MUX Voltage Ch.1 [V]'='MUX Voltage Ch.1 [V]' ORDER BY time DESC LIMIT 1 ")
        timestamp = MUX_1["TitechCoolingBox02"]["MUX Voltage Ch.1 [V]"].index[0].to_pydatetime()
        MUX_1 = MUX_1["TitechCoolingBox02"]["MUX Voltage Ch.1 [V]"][0]
        if timestamp > configChangeDate + datetime.timedelta(seconds=5):
            flag = False
        else:
            time.sleep(1)
    
    flag = True
    while flag:
        MUX_2 = DFclient.query("SELECT * FROM TitechCoolingBox02 WHERE 'MUX Voltage Ch.2 [V]'='MUX Voltage Ch.2 [V]' ORDER BY time DESC LIMIT 1 ")
        timestamp = MUX_2["TitechCoolingBox02"]["MUX Voltage Ch.2 [V]"].index[0].to_pydatetime()
        MUX_2 = MUX_2["TitechCoolingBox02"]["MUX Voltage Ch.2 [V]"][0]
        if timestamp > configChangeDate + datetime.timedelta(seconds=5):
            flag = False
        else:
            time.sleep(1)
    
    flag = True
    while flag:
        MUX_3 = DFclient.query("SELECT * FROM TitechCoolingBox02 WHERE 'MUX Voltage Ch.3 [V]'='MUX Voltage Ch.3 [V]' ORDER BY time DESC LIMIT 1 ")
        timestamp = MUX_3["TitechCoolingBox02"]["MUX Voltage Ch.3 [V]"].index[0].to_pydatetime()
        MUX_3 = MUX_3["TitechCoolingBox02"]["MUX Voltage Ch.3 [V]"][0]
        if timestamp > configChangeDate + datetime.timedelta(seconds=5):
            flag = False
        else:
            time.sleep(1)
    
    flag = True
    while flag:  
        MUX_4 = DFclient.query("SELECT * FROM TitechCoolingBox02 WHERE 'MUX Voltage Ch.4 [V]'='MUX Voltage Ch.4 [V]' ORDER BY time DESC LIMIT 1 ")
        timestamp = MUX_4["TitechCoolingBox02"]["MUX Voltage Ch.4 [V]"].index[0].to_pydatetime()
        MUX_4 = MUX_4["TitechCoolingBox02"]["MUX Voltage Ch.4 [V]"][0]
        if timestamp > configChangeDate + datetime.timedelta(seconds=5):
            flag = False
        else:
            time.sleep(1)
    return [MUX_1,MUX_2,MUX_3,MUX_4]

# enable Monitor MUX
subprocess.run(["bin/write-register","-r",controller,"-c",connectivity,"MonitorEnable","1"])
upload_conf(1,None)

for value in MUX_list:
    subprocess.run(["bin/write-register","-r",controller,"-c",connectivity,"MonitorV",str(value)])
    upload_conf(None,value)
    configChangeDate = datetime.datetime.now(datetime.timezone.utc) 
    time.sleep(1)
    voltage=getMUX(configChangeDate)
    print(voltage)
    result["Chip1"][str(value)]=voltage[0]
    result["Chip2"][str(value)]=voltage[1]
    result["Chip3"][str(value)]=voltage[2]
    result["Chip4"][str(value)]=voltage[3]

print(result)
with open("MUX result.json","w") as file:
    json.dump(result,file,indent=4)

# finish
subprocess.run(["bin/write-register","-r",controller,"-c",connectivity,"MonitorEnable","0"])
upload_conf(0,None)
DFclient.close()
client.close()

