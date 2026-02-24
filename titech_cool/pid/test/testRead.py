#!/usr/bin/env python3

import time
import datetime
import json
import traceback
import pprint
import os

#--------------------------------------------------
# Influxdb client
#--------------------------------------------------
from influxdb import InfluxDBClient
from ArduinoSingleBoxDcsController import ArduinoSingleBoxDcsController
from PeltierControl import TexioPFRControl

pid_path=os.environ["PIDPATH"]
with open(f"{pid_path}/config/default.json") as f:
    default_config = json.load(f)

cfg_file = default_config.get("config")
with open(f"{pid_path}/{cfg_file}") as f:
    config = json.load(f)


def recordToDB( ctrl = None, peltCtrl = None, host = "192.168.10.2", port = 8086, database = "REPIC", measurement = 'CoolingBox03' ):
       
    try:
        client = InfluxDBClient( host = host, port = port, database = database )
        
        if peltCtrl != None:
            data_ = [ {'fields': { 'Temperature (ch.Module) [C]'       : ctrl.m_temps['module'],
                                   'Temperature (ch.Chuck) [C]'        : ctrl.m_temps['head'],
                                   'Temperature (ch.Case) [C]'         : ctrl.m_temps['case'],
                                   'Temperature (ch.Carrier in) [C]'   : ctrl.m_temps['carrier'],
                                   'Temperature (ch.Chiller) [C]'      : ctrl.m_temps['chiller'],
                                   'Temperature (ch.Sink) [C]'         : ctrl.m_temps['sink'],
                                   'Humidity (ch.Carrier in) [%H]'     : ctrl.m_rh['carrier'],
                                   'Dew Point (ch.Carrier in) [C]'     : ctrl.m_dp['carrier'],
                                   'Hardware Interlock Status'         : ctrl.m_interlockStatus,
                                   'Interlock Status'                  : ctrl.m_interlockStatus,
                                   'Switch Heater'                     : ctrl.m_rState['heater'],
                                   'Switch LV'                         : ctrl.m_rState['lowVoltage'],
                                   'Switch HV'                         : ctrl.m_rState['highVoltage'],
                                   'Switch Chiller'                    : ctrl.m_rState['chiller'],
                                   'Switch Peltier'                    : ctrl.m_rState['peltier'],
                                   'Switch Lock'                       : ctrl.m_rState['coolingBoxLock'],
                                   'Switch Lock1'                      : ctrl.m_rState['lockState1'],
                                   'Switch Lock2'                      : ctrl.m_rState['lockState2'],
                                   'isPidActivated (ch.sum)'           : 0,
                                   'Chiller Alert'                     : ctrl.m_rState['chillerAlert'],
                                   'isInterlockEnabled (ch.sum)'       : int( ctrl.m_rState['lockState2'] == 1 ),
                                   'isInterlocked (ch.sum)'            : 0,
                                   'isInterlocked HighTemp(ch.sum)'    : int( ctrl.m_interlockStatus in [1, 2, 3] ),
                                   'isInterlocked HighDP(ch.sum)'      : int( ctrl.m_interlockStatus in [4, 5] ),
                                   'Set. Voltage (ch.Peltier) [V]'     : (1 if ctrl.m_rState['pelPlus'] == 0 else -1 )* peltCtrl.getV()['set'],
                                   'Meas. Voltage ReadBack (ch.Peltier) [V]'     : (1 if ctrl.m_rState['pelPlus'] == 0 else -1 )* peltCtrl.getV()['meas'],
                                   'Current (ch.Peltier) [A]'          : (1 if ctrl.m_rState['pelPlus'] == 0 else -1 )* peltCtrl.getI()['meas']
                               },
                   'measurement': measurement } ]
        else:
            data_ = [ {'fields': { 'Temperature (ch.Module) [C]'       : ctrl.m_temps['module'],
                                   'Temperature (ch.Chuck) [C]'        : ctrl.m_temps['head'],
                                   'Temperature (ch.Case) [C]'         : ctrl.m_temps['case'],
                                   'Temperature (ch.Carrier in) [C]'   : ctrl.m_temps['carrier'],
                                   'Temperature (ch.Chiller) [C]'      : ctrl.m_temps['chiller'],
                                   'Temperature (ch.Sink) [C]'         : ctrl.m_temps['sink'],
                                   'Humidity (ch.Carrier in) [%H]'     : ctrl.m_rh['carrier'],
                                   'Dew Point (ch.Carrier in) [C]'     : ctrl.m_dp['carrier'],
                                   'Hardware Interlock Status'         : ctrl.m_interlockStatus,
                                   'Interlock Status'                  : ctrl.m_interlockStatus,
                                   'Switch Heater'                     : ctrl.m_rState['heater'],
                                   'Switch LV'                         : ctrl.m_rState['lowVoltage'],
                                   'Switch HV'                         : ctrl.m_rState['highVoltage'],
                                   'Switch Chiller'                    : ctrl.m_rState['chiller'],
                                   'Switch Peltier'                    : ctrl.m_rState['peltier'],
                                   'Switch Lock'                       : ctrl.m_rState['coolingBoxLock'],
                                   'Switch Lock1'                      : ctrl.m_rState['lockState1'],
                                   'Switch Lock2'                      : ctrl.m_rState['lockState2'],
                                   'isPidActivated (ch.sum)'           : 0,
                                   'isInterlockEnabled (ch.sum)'       : int( ctrl.m_rState['lockState2'] == 1 ),
                                   'isInterlocked (ch.sum)'            : int( ctrl.m_interlockStatus != 0 ),
                                   'isInterlocked HighTemp(ch.sum)'    : int( ctrl.m_interlockStatus in [1, 2, 3] ),
                                   'isInterlocked HighDP(ch.sum)'      : int( ctrl.m_interlockStatus in [4, 5] )
                               },
                   'measurement': measurement } ]
        
        pprint.pprint( data_ )
        now = datetime.datetime.now()
        print( f'{now}: recorded measurements' )
        
        res = client.write_points(data_)
        
    except:
        print( traceback.format_exc() )



if __name__ == '__main__':

    while True:
        
        ctrl = ArduinoSingleBoxDcsController( { "port" : "/dev/ttyUSB_ArduinoSingleBoxController", "baud" : 9600 } )
        
        try:
            peltCtrl = TexioPFRControl( { "port":"/dev/ttyUSB_TexioPFR", "address":1, "baud":9600 } )
        except:
            peltCtrl = None
            
        
        try:
            ctrl.readAll()

            recordToDB( ctrl,
                        peltCtrl,
                        host = config.get("influxDB").get("hostname"),
                        port = config.get("influxDB").get("port"),
                        database = config.get("influxDB").get("database"),
                        measurement = config.get("id") )
            
        except:
            ctrl.reconnect()
            
        del ctrl
        del peltCtrl
        
        time.sleep(1)
