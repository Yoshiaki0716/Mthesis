#/usr/bin/python3

import time
import datetime
import matplotlib.pyplot as plt
import numpy as np
import math
import json
import traceback
import pprint

#--------------------------------------------------
# Influxdb client
#--------------------------------------------------
from influxdb import InfluxDBClient


class PidMonitor:
    #--------------------------------------------------
    def __init__(self, size):
        self.t0 = float( time.time() )
        self.size  = size
        self.data  = {}
        self.overflowStatus = 'normal'
        
    #--------------------------------------------------
    def append( self, d ):

        for key, value in d.items():
            self.data[key].append( value )

            while len(self.data[key]) > self.size:
                self.data[key].pop( 0 )
    
    #--------------------------------------------------
    def getElapsed(self):
        return float( time.time() ) - self.t0

    #--------------------------------------------------
    def printOut(self, last=-1):
        if last >=len( self.data['elapsed'] ) or last < 0:
            last = len( self.data['elapsed'] )-1
        
        form = '  |  '.join( ['{}',
                              'T module: {:6.2f}, chuck: {:6.2f}, case: {:6.2f}, carrier: {:6.2f}, DP(carrier): {:6.2f}, '
                              'P = {:7.4f}, I = {:7.4f}, D = {:7.4f}',
                              'dV = {:7.4f}',
                              'Vout = {:8.4f}  {}  {}' ] )
        
        outstr = lambda form, data, overflowStatus : form.format(
            datetime.datetime.fromtimestamp( self.t0 + data['elapsed'] ).strftime('%Y-%m-%d-%H:%M:%S'),
            data['T_module'] if data['T_module'] > -60 else None,
            data['T_chuck'],
            data['T_case'],
            data['T_carrier'],
            data['DP_carrier'],
            data['Vp'],
            data['Vi'],
            data['Vd'],
            data['dV'],
            data['V'],
            '**' + overflowStatus + '**' if overflowStatus!='normal' else '',
            'STABLE' if data['isStable'] == 1 else ''
        )

        if len(self.data) == 0:
            print( 'data is still empty' )
            return

        for n in range(-last, 0):
            try:
                d = {}
                for key in self.data.keys():
                    d.update( { key:self.data[key][n] } )
                    
                print( outstr( form, d, self.overflowStatus ) )
                
            except:
                pass
        

    #--------------------------------------------------
    def toJson(self, last_elapsed):
        if not 'elapsed' in self.data:
            return json.dumps( [] )

        i = 0
        if last_elapsed == 'all':
            i = max( 0, len(self.data)-100)
        else:
            while i < len( self.data['elapsed'] ):
                if self.data['elapsed'][i] < math.ceil( last_elapsed ):
                    i += 1
                else:
                    break
            
        dataset = { key: value[i:] for (key, value) in self.data.items() }
        
        return json.dumps( [ self.t0, dataset, self.overflowStatus ] )

    #--------------------------------------------------
    def writeLog(self, logOutFile):
        with open( logOutFile, "a" ) as F:
            dt = datetime.datetime.fromtimestamp( self.t0 + self.data['elapsed'][-1] )
            s = " ".join( [ dt.strftime('%Y-%m-%d-%H:%M:%S'),
                            str(self.t0 + self.data['elapsed'          ][-1]),
                            "{:.3f}".format(self.data['T_module'       ][-1] if self.data['T_module'       ][-1] > -60 else None ),
                            "{:.3f}".format(self.data['T_chuck'        ][-1]),
                            "{:.3f}".format(self.data['T_case'         ][-1]),
                            "{:.3f}".format(self.data['T_carrier'      ][-1]),
                            "{:.3f}".format(self.data['RH_carrier'     ][-1]),
                            "{:.3f}".format(self.data['DP_carrier'     ][-1]),
                            "{:.3f}".format(self.data['V'              ][-1]),
                            "{:.3f}".format(self.data['Vp'             ][-1]),
                            "{:.3f}".format(self.data['Vi'             ][-1]),
                            "{:.3f}".format(self.data['Vd'             ][-1]),
                            str(self.data['isStable'                   ][-1]),
                            str(self.data['isInterlockEnabled'         ][-1]),
                            str(self.data['isInterlocked'              ][-1]),
                            str(self.data['isPidActivated'             ][-1])  ] )
            F.write( s + "\n" )
            
    #--------------------------------------------------
    def recordToDB(self, logger, host = "192.168.10.7", port = 8086, database = "dcsDB", measurement = 'CpuFanCoolingBox' ):
        
        try:
            client = InfluxDBClient( host = host, port = port, database = database )
            
            data_ = [ {'fields': { 'Temperature (ch.target) [C]'       : self.data['T_set_regulated'       ][-1],
                                   'Temperature (ch.Module) [C]'       : self.data['T_module'              ][-1] if self.data['T_module'              ][-1] > -60 else None,
                                   'Temperature (ch.Chuck) [C]'        : self.data['T_chuck'               ][-1],
                                   'Temperature (ch.Case) [C]'         : self.data['T_case'                ][-1],
                                   'Temperature (ch.Carrier in) [C]'   : self.data['T_carrier'             ][-1],
                                   'Temperature (ch.Chiller) [C]'      : self.data['T_chiller'             ][-1],
                                   'Temperature (ch.Sink) [C]'         : self.data['T_sink'                ][-1],
                                   'Humidity (ch.Carrier in) [%H]'     : self.data['RH_carrier'            ][-1],
                                   'Dew Point (ch.Carrier in) [C]'     : self.data['DP_carrier'            ][-1],
                                   'Interlock Lower Limit [C]'         : self.data['DP_carrier'            ][-1]+5.0,
                                   'Interlock Upper Limit [C]'         : 35.0,
                                   'Hardware Interlock Status'         : self.data['interlockStatus'       ][-1],
                                   'Switch Heater'                     : self.data['RS_heater'             ][-1],
                                   'Switch LV'                         : self.data['RS_lowVoltage'         ][-1],
                                   'Switch HV'                         : self.data['RS_highVoltage'        ][-1],
                                   'Switch Chiller'                    : self.data['RS_chiller'            ][-1],
                                   'Switch Peltier'                    : self.data['RS_peltier'            ][-1],
                                   'Switch Lock'                       : self.data['RS_coolingBoxLock'     ][-1],
                                   'Switch Lock1'                      : self.data['RS_lockState1'         ][-1],
                                   'Switch Lock2'                      : self.data['RS_lockState2'         ][-1],
                                   'Chiller Alert'                     : self.data['RS_chillerAlert'       ][-1],
                                   'Set. Voltage (ch.Peltier) [V]'     : self.data['V'                     ][-1],
                                   'Set. Voltage ReadBack (ch.Peltier) [V]'     : self.data['setV'         ][-1],
                                   'Meas. Voltage ReadBack (ch.Peltier) [V]'    : self.data['measV'        ][-1],
                                   'Current (ch.Peltier) [A]'          : self.data['I'                     ][-1],
                                   'PID Feedback (ch.P) [V]'           : self.data['Vp'                    ][-1],
                                   'PID Feedback (ch.I) [V]'           : self.data['Vi'                    ][-1],
                                   'PID Feedback (ch.D) [V]'           : self.data['Vd'                    ][-1],
                                   'PID Feedback (ch.sum) [V]'         : self.data['dV'                    ][-1],
                                   'TemperatureSlope'                  : float( self.data['slope'                 ][-1] ),
                                   'isConverging (ch.sum)'             : self.data['isConverging'          ][-1],
                                   'isStable (ch.sum)'                 : self.data['isStable'              ][-1],
                                   'isPidActivated (ch.sum)'           : self.data['isPidActivated'        ][-1],
                                   'isInterlockEnabled (ch.sum)'       : self.data['isInterlockEnabled'    ][-1],
                                   'isInterlocked (ch.sum)'            : self.data['isInterlocked'         ][-1],
                                   'isInterlocked HighTemp(ch.sum)'    : self.data['isInterlockedHighTemp' ][-1],
                                   'isInterlocked HighDP(ch.sum)'      : self.data['isInterlockedHighDP'   ][-1],
                                   'isInterlocked IrregularTargetTemp(ch.sum)'      : self.data['isInterlockedIrregularTargetTemp'   ][-1],
                                   'isInterlocked PeltierTripped(ch.sum)'      : self.data['isInterlockedPeltierTripped'   ][-1],
                                   'Interlock Status'                  : self.data['interlockStatus'       ][-1],
                                   'Temperature RMS (ch.target) [K]'   : self.data['T_rms'                 ][-1],
                                   'Voltage RMS (ch.target) [V]'       : self.data['V_rms'                 ][-1]   },
                       'measurement': measurement } ]
            
            # pprint.pprint( data_ )
            
            res = client.write_points(data_)
            
        except:
            logger.error( pprint.pformat( traceback.format_exc() ) )
            logger.error('PidMonitor.recordToDB(): failure in DB recording!')
            

            
