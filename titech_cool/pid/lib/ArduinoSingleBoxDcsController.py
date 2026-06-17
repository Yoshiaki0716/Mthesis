#/usr/bin/python3

import serial
import time
import os
import numpy
import pprint
import datetime
import statistics
import copy
from multiprocessing import Process, Value

class ArduinoSingleBoxDcsController:
    m_port = '/dev/ttyUSB_ArduinoSingleBoxController'
    m_baud = 9600
    m_temps = { 'module': -9999., 'chiller': -9999., 'case' : -9999., 'sink' : -9999., 'head' : -9999., 'carrier' : -9999. }
    m_temps_raw = { 'chiller': -9999., 'case' : -9999., 'sink' : -9999., 'head' : -9999. }
    m_temps_ref = { 'chiller': -9999., 'case' : -9999., 'sink' : -9999., 'head' : -9999. }
    m_rh = { 'carrier' : -9999. }
    m_dp = { 'carrier' : -9999. }
    m_rState = {'heater' : -1, 'peltier' : -1, 'chiller' : -1, 'pelPlus' : -1, 'pelMinus' : -1,
                'lowVoltage' : -1, 'highVoltage' : -1, 'coolingBoxLock' : -1, 'lockState1': -1, 'lockState2': -1, 'chillerAlert': -1 }
    m_interlockStatus = -1
    m_busy = False

    m_useCalib = True


    '''
    m_calibMap = { 'chiller' : [ 3.04792, 1.04162, -0.0052014 ],
                   'sink'    : [ 3.69336, 1.07189, -0.00558717 ],
                   'case'    : [ 3.04063, 1.09587, -0.00383257 ],
                   'head'    : [ 2.44063, 1.0517, -0.00229164 ] }
    
    # adjustment at KEK
    m_calibMap = { 'chiller' : [ 2.00, 1.04162, -0.0052014 ],
                   'sink'    : [ 1.60, 1.07189, -0.00558717 ],
                   'case'    : [ -0.25, 1.09587, -0.00383257 ],
                   'head'    : [ 0.20, 1.0517, -0.00229164 ] }
    '''
    
    # re-calibrated: Feb 8 2022
    m_calibMap = { 'chiller' : [ 2.00, 1.04162, -0.0052014 ],
                   'sink'    : [ 1.60, 1.07189, -0.00558717 ],
                   'case'    : [ 0.60, 1.09587, -0.00383257 ],
                   'head'    : [ -0.267, 1.074, -0.002736 ] }

    
    def __init__(self, config):
        if 'port' in config:
            self.m_port = config['port']
        if 'baud' in config:
            self.m_baud = config['baud']
        if 'useCalib' in config:
            self.m_useCalib = config['useCalib']

        os.system('stty -F ' + self.m_port + ' -hupcl')
        
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        time.sleep(0.3)

        self.m_busy = False
        self.m_device_error = False

    def is_device_ok(self):
        return self.m_device_error

    def calib( self, temp, temp_ref, channel ):
        if self.m_useCalib:
            calib = self.m_calibMap[channel][0] + self.m_calibMap[channel][1] * temp + self.m_calibMap[channel][2] * temp * temp
            return calib

        else:
            return temp

    def waitBusy(self):
        while self.m_busy:
            time.sleep(0.01)
            
    def reconnect(self):
        self.m_ctrl.close()
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        time.sleep(1)
        
        
    def command( self, string ):
        
        while True:
            self.m_busy = True
            
            self.m_ctrl.write( (string+'\n').encode() )
            
            try:
                out = self.m_ctrl.readline().decode().rstrip()
            except Exception as e:
                print( f'WARNING: command( {string} ) failed: exception {str(e)}' )
                self.m_busy = False
                self.m_device_error = True

            if out.find('Invalid') == 0:
                
                self.m_ctrl.write( '\n'.encode() )
                self.m_ctrl.readline().decode().rstrip()
                
                continue
                
            self.m_busy = False
            
            break
                
        return out
        

    def temps(self):
        
        self.waitBusy()
        
        
        # keep the previous value
        temps_prev = copy.deepcopy( self.m_temps )
        
        out = self.command( 'Temps' ).split(',')
        
        try:
            ts = [ float(val) for val in out ]
        except Exception as e:
            print( f'WARNING: temps() failed: exception {str(e)}' )
            self.m_busy = False
            return
        
        while len( ts ) < 8 :
            out = self.m_ctrl.readline().decode().rstrip().split(',')
            ts = [ float(val) for val in out ]
            
        # here, temp channels are defined as ['chiller', 'case', 'sink', 'head' ]
    
        self.m_temps['chiller'] = self.calib( ts[0], ts[4], 'chiller' )
        self.m_temps['case']    = self.calib( ts[2], ts[6], 'case' )
        self.m_temps['sink']    = self.calib( ts[1], ts[5], 'sink' )
        self.m_temps['head']    = self.calib( ts[3], ts[7], 'head' )
        
        self.m_temps_raw['chiller'] = ts[0]
        self.m_temps_raw['case']    = ts[2]
        self.m_temps_raw['sink']    = ts[1]
        self.m_temps_raw['head']    = ts[3]
        
        self.m_temps_ref['chiller'] = ts[4]
        self.m_temps_ref['case']    = ts[6]
        self.m_temps_ref['sink']    = ts[5]
        self.m_temps_ref['head']    = ts[7]
        
        for k,v in self.m_temps.items():
            if k == 'module': continue
            if numpy.isnan(v):
                self.m_temps[k] = -273.15
    
        temp_array = []
        for i in range(9):
            temp = self.m_temps['module']
            out = self.command( 'ModuleTemp' )
            self.m_busy = False
            
            if out == 'IRREGULAR':
                self.m_temps['module'] = None
            else:
                try:
                    temp = float( out )
                except:
                    self.m_temps['module'] = None
                    
            if temp is not None and temp < -60:
                temp = None

            if temp is not None:
                temp_array.append( temp )
            
        self.m_temps['module'] = sum( temp_array ) / (1.0 * len(temp_array) ) if len( temp_array ) == 9 else None
        
        
    def carrierEnv(self):
        self.waitBusy()
        
        for i in range(3):
            
            try:
                out = self.command( 'CarrierEnv' ).split(',')
                vals = [ float(val) for val in out ]
            except Exception as e:
                print( f'WARNING: carrierEnv() failed: exception {str(e)}' )
                self.m_busy = False
                return
                
                # here, vals[0] = temperature, vals[1] = relative humidity
            if ( abs( self.m_temps['carrier'] - vals[0] ) > 1.0 or
                 abs( self.m_rh['carrier'] - vals[1] ) > 1.0 or
                 abs( self.m_dp['carrier'] - vals[2] ) > 1.0  ):
                continue

            else:
                break
            
        self.m_temps['carrier'] = vals[0]
        self.m_rh['carrier']    = vals[1]
        self.m_dp['carrier']    = vals[2]
        
    def relayOn( self, channel ):
        if( channel < 1 or channel > 8 ):
            raise ValueError('Invalid channel number!')

        self.waitBusy()
        
        try:
            out = self.command( 'RelayOn ' + str(channel) )
        except Exception as e:
            print( f'WARNING: relayOff() failed: exception {str(e)}' )
            self.m_busy = False
            return
        
    def relayOff( self, channel ):
        if( channel < 1 or channel > 8 ):
            raise ValueError('Invalid channel number!')
        
        self.waitBusy()
        
        try:
            out = self.command( 'RelayOff ' + str(channel) )
        except Exception as e:
            print( f'WARNING: relayOff() failed: exception {str(e)}' )
            self.m_busy = False
            return
        
    def relayStates(self):

        self.waitBusy()
        
        try:
            out = self.command( 'RelayStates' ).split(',')
        except Exception as e:
            print( f'WARNING: relayStates() failed: exception {str(e)}' )
            self.m_busy = False
            return
        
        try:
            vals = [ 1 if int(val)>0 else 0 for val in out ]
        except Exception as e:
            print( f'WARNING: relayStates() failed: exception {str(e)}' )
            self.m_busy = False
            return
        
        self.m_rState['heater']         = vals[0]
        self.m_rState['peltier']        = vals[1]
        self.m_rState['chiller']        = vals[2]
        self.m_rState['pelPlus']        = vals[3]
        self.m_rState['pelMinus']       = vals[4]
        self.m_rState['lowVoltage']     = vals[5]
        self.m_rState['highVoltage']    = vals[6]
        self.m_rState['coolingBoxLock'] = vals[7]
        self.m_rState['lockState1']     = vals[8]
        self.m_rState['lockState2']     = vals[9]
        self.m_rState['chillerAlert']   = vals[10]
            
        
    def pelPolNormal(self):
        self.waitBusy()
        
        out = self.command( 'PelPolNormal' )

        self.relayStates()
        
    def pelPolInverted(self):
        self.waitBusy()
        
        out = self.command( 'PelPolInverted' )

        self.relayStates()

    
    def interlockStatus(self):
        self.waitBusy()
        
        out = self.command( 'Status' )
        
        try:
            self.m_interlockStatus =  int(out) 
        except Exception as e:
            self.m_busy = False
            print( f'WARNING: interLockStatus() failed: exception {str(e)}' )
            return

################################################################################################
    def interlockTest(self,mode):
        if(mode<0 or mode>7):
            raise ValueError('Invalid mode number')
        self.waitBusy()

        try:
            out = self.command('InterlockTest' + str(mode))
        except Exception as e:
            print(f'WARNING: interlockTest() failed: exception {str(e)}')
            self.m_busy = False
            return
    ################################################################################

    def unlock(self):
        self.waitBusy()
        
        out = self.command( 'Unlock' )
        self.relayStates()

    def reset(self):
        self.waitBusy()
        
        out = self.command( 'Reset' )
        self.interlockStatus()
        self.relayStates()
        
    def lockStatus(self):
        self.relayStates()
        #return ( self.m_rState['lockState1'] == 0 or self.m_rState['lockState2'] == 0 ) # temporary
        return self.m_rState['lockState2'] == 0 # temporary
        #return ( self.m_rState['lockState1'] == 0 and self.m_rState['lockState2'] == 0 )
        
    def readAll(self):
        self.temps()
        self.carrierEnv()
        self.relayStates()
        self.interlockStatus()

    def print(self):
        interlockModes = { 0 : 'Normal', 1 : 'HighTemp, LV/HV down', 2 : 'HighTemp, LV/HV and Peltier down',
                           3 : 'HighTemp, Idling', 4 : 'HighDewPoint, Peltier/Chiller down', 5 : 'HighDewPoint, Idling', 6 : 'Chiller Alert', 7 : 'NTC Disconnected' }
        print( 'Interlock: ', interlockModes[self.m_interlockStatus] )
        print( 'Relays: ' )
        pprint.pprint( self.m_rState )
        print( 'Temps: ' )
        pprint.pprint( self.m_temps )
        print( 'Temps_raw: ' )
        pprint.pprint( self.m_temps_raw )
        print( 'Temps_Ref: ' )
        pprint.pprint( self.m_temps_ref )
        print( 'RH: ', self.m_rh )
        print( 'DP: ', self.m_dp )

prevTime = datetime.datetime.now()
totalTime = prevTime - prevTime

def trackTime( title='', reset = False ):
    
    global prevTime, totalTime
    
    now = datetime.datetime.now()
    
    if not reset:
        totalTime += (now - prevTime)
    else:
        totalTime -= totalTime
    
    print( f'{title}: lap = {(now - prevTime).total_seconds()*1000:.1f} ms, total = {totalTime.total_seconds()*1000:.1f} ms' )
    prevTime = now


if __name__ == '__main__':
    
    trackTime('begin')
    ctrl = ArduinoSingleBoxDcsController( { "port" : "/dev/ttyUSB_ArduinoSingleBoxController", "baud" : 9600, 'useCalib':True } )
    
    out = ctrl.m_ctrl.readline().decode().rstrip()
    print(f"residual buffer: '{out}'")
    time.sleep(1)
    
    ctrl.temps()
    ctrl.carrierEnv()
    ctrl.relayOn(2)
    ctrl.relayOff(2)
    ctrl.relayStates()
    ctrl.pelPolInverted()
    ctrl.pelPolNormal()
    ctrl.interlockStatus()
    #ctrl.reset()
    ctrl.print()
    
    print('\n\nTime challenge:')
    meas = []
    
    for i in range(10):
        print( '----------------------' )
        trackTime( 'init', True )
        ctrl.readAll()
        trackTime( 'readAll' )
        meas += [ totalTime.total_seconds()*1000 ]
        ctrl.print()
        #trackTime( 'print', True )
    
    print( 'Average: {:.1f} ms, Stdev: {:.1f} ms'.format( statistics.mean( meas ), statistics.stdev( meas ) ) )
