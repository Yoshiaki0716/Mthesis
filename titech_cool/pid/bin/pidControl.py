#!/usr/bin/python

import sys
from sys import argv

import json
from collections import OrderedDict 

from influxdb import InfluxDBClient

import RPi.GPIO as GPIO
import threading

from statistics import mean

from PeltierControl import TakasagoKXControl
from PeltierControl import TakasagoZXControl

from MAX31855 import *
from MCP3208  import *
from ArduinoSHT85Control import *

from PidControl import *
from PidMonitor import *


#--------------------------------------------------------------------------------




#--------------------------------------------------------------------------------
def configure():
    global pid
    global peltCtrl
    
    argvs = sys.argv
    argc = len(argvs)
    
    if argc < 3:
        print("Something is wrong. Usage: python Control.py Takasago.json logOutFile")
        sys.exit()
        
    
    jsonfilename = sys.argv[1]
    pid.logOutFile   = 'logs/' + sys.argv[2]
    
    #--------------------------------------------------
    # reading configs
    #--------------------------------------------------
    
    with open(jsonfilename) as f:
        j = json.load(f, object_pairs_hook=OrderedDict)
        
        pid.T_target    = float( j["T_target"] )
        pid.K_p         = j["Pcontrol"]
        pid.K_i         = j["Icontrol"]
        pid.K_d         = j["Dcontrol"]
        pid.interval    = j["interval"]
        pid.historySize = j["historySize"]
        pid.minV        = j["minV"]
        pid.maxV        = j["maxV"]
        
        if j["devicetype"] == "TakasagoKX":
            peltCtrl = TakasagoKXControl( j["port"], j["address"] )
        elif j["devicetype"] == "TakasagoZX":
            peltCtrl = TakasagoZXControl( j["port"], j["address"] )
    
    
    if argc >= 4:
        if( float( sys.argv[3] ) > -99 ):
            peltCtrl.setV( float( sys.argv[3] ) )
        
    if argc >= 5:
        pid.T_target = float( sys.argv[4] )
    
    
#--------------------------------------------------------------------------------
if __name__ == '__main__':

    #--------------------------------------------------
    # instantiations
    #--------------------------------------------------
    
    pid = PidControl()
    mon = PidMonitor(500)

    #sht85 = ArduinoSHT85Control("/dev/ttyUSB2")
    #sht85 = ArduinoSHT85Control("/dev/ttyACM2") 
    sht85 = ArduinoSHT85Control( { "port" : "/dev/ttyACM2", "baud" : 9600 } ) 
    
    thermistor = ThermistorReader( device = MAX31855(11),
                                   channel = 0,
                                   calib = ThermistorCalibConfig( [ 0.0006892, -0.05584, 0.99 ] ) )
    
    ntc = NtcReader( device  = MCP3208( ss = 7, clk = 33, miso = 35, mosi = 37, w = 1.e-3),
                     channel = 0,
                     calib   = NtcCalibConfig( R25 = 100, Rext = 470, B = 4281, T0 = 273.15+27.16-1.5) )
    
    sensors  = { 'base' : thermistor, 'surface' : ntc, 'air' : sht85 }
    peltCtrl = None
    
    #--------------------------------------------------
    # configurations
    #--------------------------------------------------
    
    configure()
    
    pid.injectDependencies( mon, sensors, peltCtrl )
    
    peltCtrl.setOn()
    pid.V = peltCtrl.getV()
    
    
    #--------------------------------------------------
    # start process threads
    #--------------------------------------------------
    
    tPID = threading.Thread( target = pid.pidLoop )
    tPID.start()
    
    tPlot = threading.Thread( target = mon.plot, args = ( [ pid ] ) )
    tPlot.start()
    
    tCommand = threading.Thread( target = pid.userCommand )
    tCommand.start()
    
    
    #--------------------------------------------------
    # end process
    #--------------------------------------------------
    
    tCommand.join()
    tPlot.join()
    tPID.join()
    
    peltCtrl.setV(0)
    peltCtrl.setOff()

'''
    #influxdb
    #client = InfluxDBClient(host="192.168.1.44",port=8086,database="testdb")
    #for KEK (to obsidian)
    client = InfluxDBClient(host="192.168.10.7",port=8086,database="testdb")                                                     
    #data_ = [{'fields': {'Termo1': T_measure[1]},'measurement': 'testms','tags': {'cat1': 'titech'}}]
    data_ = [{'fields': {'Thermo': Ts[N-1] , 'V_Peltier': V , 'P': P , 'I' : I , 'D' : D },'measurement': 'testms','tags': {'cat1': 'titech'}}]
    res = client.write_points(data_)
'''
