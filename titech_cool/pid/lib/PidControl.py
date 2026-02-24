import os
import sys
import socket
import time
import math
import numpy as np
from statistics import mean
from functools import reduce
from functools import partial
from operator import add
import json
import traceback
import pprint
import multiprocessing
from Kalman import KalmanSmoother
from DeviceInvoker import *

import smtplib
import datetime
from email.mime.text import MIMEText
import subprocess
import threading

from skpy import Skype
from skpy.msg import SkypeMsg

from PeltierControl import TexioPFRControl

#--------------------------------------------------------------------------------
def getSlope( dTs, ts ):
    span = min(max(10, 30 - int(abs(dTs[-1]))), len(ts))

    x = np.array( [ t - ts[-span] for t in ts[-span:] ] )
    y = np.array( KalmanSmoother( dTs[-span:], 0.05, 0.3 ).process() )
    weights = np.array( [ 1./0.25 for dT in dTs[-span:] ] )
    
    params = None
    
    if len(dTs) < 5:
        return 0
    else:
        params = np.polyfit( x, y, 2, w = weights )
    
    f = np.poly1d( params )
    d = np.polyder( f )

    return d( x[-1] )
    

#--------------------------------------------------------------------------------
def multiMeas( meas, n ):
    while n > 0:
        yield meas()
        n -= 1


#--------------------------------------------------------------------------------
class PidControl:
    def kill(self, args):
        self.kill = True
        self.conn.sendall( ('accepted: kill' + "\n").encode() )

    def params(self, args):
        j = json.dumps( { 'T_set'              : self.T_set,
                          'T_set_regulated'    : self.T_set_regulated,
                          'V'                  : self.V,
                          'Kp'                 : self.Kp,
                          'Ki'                 : self.Ki,
                          'Kd'                 : self.Kd,
                          'maxPfeedback'       : self.maxPfeedback,
                          'maxIfeedback'       : self.maxIfeedback,
                          'maxDfeedback'       : self.maxDfeedback,
                          'integralDuration'   : self.integralDuration,
                          'tOffsetP'           : self.tOffsetP,
                          'tOffsetI'           : self.tOffsetI,
                          'tOffsetD'           : self.tOffsetD,
                          'interval'           : self.interval,
                          'jumpTolerance'      : self.jumpTolerance,
                          'isInterlockEnabled' : self.isInterlockEnabled,
                          'isInterlocked'      : self.isInterlocked,
                          'isInterlockedHighTemp'      : self.isInterlockedHighTemp,
                          'isInterlockedHighDP'      : self.isInterlockedHighDP,
                          'isInterlockedPeltierTripped'      : self.isInterlockedPeltierTripped,
                          'isPidActivated'     : self.isPidActivated   } )
        self.conn.sendall( (str(j)+ '\n').encode()  )
    
    def execDCS( self, func ):
        
        try:
            return func()
        except Exception as e:
            self.logger.debug( pprint.pformat( e ) )
            self.logger.debug( 'Error in dcsCtrl ==> Reconnect' )
            time.sleep(1)
            self.dcsCtrl.reconnect()
            self.logger.debug( 'dcsCtrl reconnected.' )
            return func()
    
        if self.is_dcsCtrl_good and not self.dcsCtrl.is_device_ok():
            self.isPidActivated = False
            self.sendEmail( 'HARDWARE ERROR', 'Lost communication with Arduino', category = 'FATAL' )
            os.environ["DISPLAY"]=":1"
            os.system( "${EQC_OPERATOR_DIR}/lib/Notifier.py --message \"FATAL: HARDWARE_ERROR: Lost communication with Arduino\" --color red &" )
            self.logger.error('setT(): Invalid range value specified!')
            self.shutdownMPOD()
            self.is_dcsCtrl_good = False
            
        
    def setT(self, args):
        try:
            val = float( args[0] )
            
            minT    = self.config['pidControl']['minT']
            maxT    = self.config['pidControl']['maxT']

            if self.isInterlockEnabled:
                if val < minT or val > maxT:
                    self.logger.error('setT(): Invalid range value specified!')
                    self.conn.sendall( ('accepted: setT: Invalid range value specified! T must be {} < T < {}\n'.format( minT, maxT )).encode() )
                else:
                    self.T_set = float( args[0] )
                    self.conn.sendall( ('accepted: setT ' + args[0] + "\n").encode() )
                    self.sendEmail( 'SetT {}'.format( args[0] ) )
            else:
                if val < 15 or val > maxT:
                    self.conn.sendall( ('accepted: setT: Invalid range value specified! T must be {} < T < {}\n'.format( 15, maxT )).encode() )
                else:
                    self.T_set = float( args[0] )
                    self.conn.sendall( ('accepted: setT ' + args[0] + "\n").encode() )
            
        except:
            self.conn.sendall( f'setT: exception raised for setT {args[0]}\n'.encode() )

    def resetV(self, args):
        try:
            V_new = float( args[0] )
        except:
            self.conn.sendall( 'setV: invalid phrase (typo??)\n'.encode() )
            return

        try:            
            if self.V >= 0.0 and V_new < 0.0:
                self.peltCtrl.setV( 0.0 )
                time.sleep(1)
                self.execDCS( self.dcsCtrl.pelPolInverted )
                time.sleep(0.5)
                
            elif self.V < 0.0 and V_new >= 0.0:
                self.peltCtrl.setV( 0.0 )
                time.sleep(1)
                self.execDCS( self.dcsCtrl.pelPolNormal )
                time.sleep(0.5)
                
            self.V = V_new
            self.peltCtrl.setV( V_new )
            
            self.logger.warning( 'resetV: V = {} [V]'.format( self.V ) )
            self.conn.sendall( ('accepted: resetV ' + args[0] + "\n").encode() )
            
        except Exception as e:
            self.logger.error( str(e) )
            self.logger.error( pprint.pformat(traceback.format_exc()) )
            self.logger.error('Error in resetV, possibly peltier control was lost')
            self.sendEmail( 'Failure in resetV()', f'Error in resetV, possibly peltier control was lost.\n```\n{str(e)}\n```\n\n```\n{ traceback.format_exc() }\n```', category = 'ERROR' )
            return

    def setRelayOn(self, args):
        try:
            channel = int( args[0] )
            
            self.logger.warning( 'setRelayOn: {}'.format( channel ) )
            
            self.execDCS( partial( self.dcsCtrl.relayOn, int( args[0] ) ) )
                
            self.conn.sendall( ('accepted: setRelayOn ' + args[0] + "\n").encode() )
        except:
            self.conn.sendall( 'setRelayOn: invalid phrase (typo??)\n'.encode() )

    def setRelayOff(self, args):
        try:
            channel = int( args[0] )
            
            self.logger.warning( 'setRelayOff: {}'.format( channel ) )
            
            self.execDCS( partial( self.dcsCtrl.relayOff, int( args[0] ) ) )
                
            self.conn.sendall( ('acceptedwienermpod_ivi: setRelayOff ' + args[0] + "\n").encode() )
        except:
            self.conn.sendall( 'setRelayOff: invalid phrase (typo??)\n'.encode() )

    def setKp(self, args):
        try:
            self.Kp = float( args[0] )
            self.conn.sendall( ('accepted: setKp ' + args[0] + "\n").encode() )
        except:
            self.conn.sendall( 'setKp: invalid phrase (typo??)\n'.encode() )
    
    def setKi(self, args):
        try:
            self.Ki = float( args[0] )
            self.conn.sendall( ('accepted: setKi ' + args[0] + "\n").encode() )
        except:
            self.conn.sendall( 'setKi: invalid phrase (typo??)\n'.encode() )
    
    def setKd(self, args):
        try:
            self.Kd = float( args[0] )
            self.conn.sendall( ('accepted: setKd ' + args[0] + "\n").encode() )
        except:
            self.conn.sendall( 'setKd: invalid phrase (typo??)\n'.encode() )

    def setMaxP(self, args):
        try:
            self.maxPfeedback = float( args[0] )
            self.conn.sendall( ('accepted: setMaxP ' + args[0] + "\n").encode() )
        except:
            self.conn.sendall( 'setMaxP: invalid phrase (typo??)\n'.encode() )

    def setMaxI(self, args):
        try:
            self.maxIfeedback = float( args[0] )
            self.conn.sendall( ('accepted: setMaxI ' + args[0] + "\n").encode() )
        except:
            self.conn.sendall( 'setMaxI: invalid phrase (typo??)\n'.encode() )

    def setMaxD(self, args):
        try:
            self.maxDfeedback = float( args[0] )
            self.conn.sendall( ('accepted: setMaxD ' + args[0] + "\n").encode() )
        except:
            self.conn.sendall( 'setMaxD: invalid phrase (typo??)\n'.encode() )

    def setIduration(self, args):
        try:
            self.integralDuration = int( args[0] )
            self.conn.sendall( ('accepted: setIduration ' + args[0] + "\n").encode() )
        except:
            self.logger.warning( 'setIduration: invalid phrase (typo??)' )

    def setTimeOffsetP(self, args):
        try:
            self.tOffsetP = int( args[0] )
            self.conn.sendall( ('accepted: setTimeOffsetP ' + args[0] + "\n").encode() )
        except:
            self.logger.warning( 'setTimeOffsetP: invalid phrase (typo??)' )

    def setTimeOffsetI(self, args):
        try:
            self.tOffsetI = int( args[0] )
            self.conn.sendall( ('accepted: setTimeOffsetI ' + args[0] + "\n").encode() )
        except:
            self.logger.warning( 'setTimeOffsetI: invalid phrase (typo??)' )

    def setTimeOffsetD(self, args):
        try:
            self.tOffsetD = int( args[0] )
            self.conn.sendall( ('accepted: setTimeOffsetD ' + args[0] + "\n").encode() )
        except:
            self.logger.warning( 'setTimeOffsetD: invalid phrase (typo??)' )

    def restoreInterlock(self, args):
        if self.isHumiditySafe():
            self.isInterlocked = False
            self.isInterlockedHighTemp = False
            self.isInterlockedHighDP = False
            self.isInterlockedPeltierTripped = False
            self.isInterlockedIrregularTargetTemp   = False
            
            # Restore LV Relay
            self.execDCS( partial( self.dcsCtrl.relayOn, 6 ) )
            
            # Restore HV Relay
            self.execDCS( partial( self.dcsCtrl.relayOn, 7 ) )
            
            # Restore Peltier PS
            self.peltCtrl.resetAlarms()
            
            time.sleep(1)
            
            status = self.peltCtrl.setV(0.0)
            self.peltCtrl.setOff()
            time.sleep(1)
            status = self.peltCtrl.setV( self.V )
            time.sleep(3)
            self.peltCtrl.setOn()
            
            self.conn.sendall( 'Interlock was restored.\n'.encode() )
            self.sendEmail( 'INTERLOCK_RESTORED', category='WARNING' )
        else:
            self.conn.sendall( 'Could not restore since humidity is still in unsafe range.\n'.encode()  )
    
    
    def stopPID(self, args = None):
        if self.isHumiditySafe():
            if not self.isInterlockEnabled:
                objName = self.config['pidControl']['object']
                if self.mon.data[self.objName][-1] is not None and self.mon.data[self.objName][-1] > self.roomDPMon.getLast() + 5.0:
                    self.isPidActivated = False
                    self.conn.sendall( ('accepted: stopped PID feedback\n').encode() )
                    self.sendEmail( 'Stopped PID control', 'Stopped PID control by user action' )
                else:
                    self.conn.sendall( ('rejected: at this moment cannot stop PID feedback since the temperature is too low!\n').encode() )
            else:
                self.isPidActivated = False
                self.conn.sendall( ('accepted: stopped PID feedback\n').encode() )
                self.sendEmail( 'Stopped PID control', 'Stopped PID control by user action' )
        else:
            self.conn.sendall( ('rejected: not humidity safe\n').encode() )
    
    def restartPID(self, args):
        if self.isHumiditySafe() and (not self.isInterlocked):
            self.isPidActivated = True
            self.conn.sendall( 'PID control was re-activated.\n'.encode()  )
            self.sendEmail( 'Restarted PID control', 'Restarted PID control by user action' )
        else:
            self.logger.warning( 'restartPID: could not restart since humidity is still in unsafe range or the interlock is still active.' )
            self.conn.sendall( 'could not restart since humidity is still in unsafe range or the interlock is still active.\n'.encode()  )
            
    def unlock( self, args ):
        
        self.checkCoverStatus()
        if self.lockStatus == 0:
            self.conn.sendall( ('rejected: cover should be already open.\n').encode() )
            return
        
        if self.isHV_On():
            self.conn.sendall( ('rejected: cannot unlock the cover at this moment. HV is still on.\n').encode() )
            return
        
        if self.isLV_On():
            self.conn.sendall( ('rejected: cannot unlock the cover at this moment. LV is still on.\n').encode() )
            return
            
        if self.mon.data[self.altName][-1] > self.roomDPMon.getLast() + 5.0:
            self.isPidActivated = True
            self.T_set = max( self.T_default, self.roomDPMon.getLast() + 5.0 )
            self.T_set_regulated = max( self.T_default, self.roomDPMon.getLast() + 5.0 )
            
            self.execDCS( self.dcsCtrl.unlock )
            self.lockStatus = 0
            
            self.conn.sendall( ('accepted: unlocked the cover with safe. Interlock is now inactive.\n').encode() )
            self.sendEmail( 'COVER_UNLOCKED', category = 'WARNING' )
        else:
            mon_data = self.mon.data[self.altName][-1]
            mon_roomDP = self.roomDPMon.getLast()
            self.conn.sendall( (f'rejected: cannot unlock the cover at this moment due to condensation risk with cover opening. T ={mon_data}, roomDP={mon_roomDP} \n').encode() )
        
    def resetDCS( self, args ):
        self.execDCS( self.dcsCtrl.reset )
        self.conn.sendall( ('accepted: sent a reset command to Arduino.\n').encode() )
        time.sleep(3)
    
        # reset HV latch
        cmd = f'snmpset -v 2c -m +WIENER-CRATE-MIB -c guru {self.mpodIP} groupsSwitch.64 i 10'
        self.logger.debug( cmd )
        subprocess.call( cmd, shell=True )
        
        
    def snmpcmd( self, cmdType, name, channel, valType, value ):
        cmd = f'snmp{cmdType} -OqvU -v 2c -m +WIENER-CRATE-MIB -c guru {self.mpodIP} {name}.{channel} {valType} {value}'
        self.logger.debug( cmd )
        try:
            subout = subprocess.run( cmd, capture_output=True, shell=True, check=True )
            return subout.stdout.decode()
        except Exception as e:
            self.logger.warning( str(e) )
            return ""
        
    def mpodcmd( self, channel, mpodCommand, args=None ):
        mpodDir='/home/admin/wienermpod_ivi'
        mpodCfg='configs/lr_powersupply.json'
        cmd = f'cd {mpodDir}; python3 libDCS/mpodOperator.py -e {mpodCfg} -c {channel} {mpodCommand}'
        
        if args != None:
            if args[0]=='cli':
                self.conn.sendall( ( f'accepted: MPOD channel {channel} --> {mpodCommand}.\n').encode() )
            
        self.logger.info( cmd )
        stdout = subprocess.run( cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT ).stdout.decode()
        # self.logger.info( stdout )
        
    #--------------------------------------------------
    def RIGOLcmd( self, ip, channel, RIGOLCommand, args=None):
        RIGOLDir='/nas/dcs/wienermpod_ivi'
        cmd = f"cd {RIGOLDir}; python3 -c \"from RIGOL_DP821 import RIGOL_DP821_PYVISA; RIGOL_DP821_PYVISA().{RIGOLCommand}\"  {ip}"
        if args != None:
            if args[0]=='cli':
                self.conn.sendall( ( f'accepted: RIGOL channel {channel} --> {RIGOLCommand}.\n').encode() )
        self.logger.info( cmd )
        stdout = subprocess.run( cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT ).stdout.decode()
        return stdout
        
    #--------------------------------------------------
    def setLV_On( self, args=None ):
        self.logger.debug( "Called setLV_On()" )
        self.checkCoverStatus()
        if self.lockStatus == 1:
            if self.LVtype == "MPOD":
                self.mpodcmd( self.mpodLVName, 'power-on', args )
            elif self.LVtype == "RIGOL":
                self.RIGOLcmd(self.LVip, self.LVch, "turnOn()", args)
            elif self.LVtype == "CUSTOM":
               # self.RIGOLcmd(self.LVip, self.LVch, "turnOn()", args)
                if args != None:
                    if args[0]=='cli':
                        self.conn.sendall( ( f'accepted: LV to set ON.\n').encode() )
                self.logger.info( "CUSTOM LV command called: {}".format(self.LVOnCmd) )
                stdout = subprocess.run( self.LVOnCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT ).stdout.decode()
                self.logger.debug( stdout )
        else:
            if args and args[0]=='cli':
                self.conn.sendall( ('rejected: cover is open!\n').encode() )
            else:
                self.logger.warning( "setLV_On rejected: cover is open." )
        
    #--------------------------------------------------
    '''
    def setLV_Off( self, args=None ):
        
        if self.isHV_On():
            if args[0]=='cli':
                self.conn.sendall( ('rejected: HV is still ON!\n').encode() )
        elif self.LVtype == "MPOD":
            self.mpodcmd( self.mpodLVName, 'power-off', args )
        elif self.LVtype == "RIGOL":
            self.RIGOLcmd(self.LVip, self.LVch, "turnOff()", args)
        else :
          self.conn.sendall( ('rejected: LV Power Off command  !\n').encode() ) 
    '''
    # -------------------------------------------------
    def setLV_Off( self, args=None ):
        self.logger.debug( "Called setLV_Off()" )
        if args[0] != 'interlock' and self.isHV_On():
            try:
                if args[0]=='cli':
                    self.conn.sendall( ('rejected: HV is still ON!\n').encode() )
            except:
                self.logger.warning( 'rejected: HV is still ON!' )

        else :
            if self.LVtype == "MPOD":
                self.mpodcmd( self.mpodLVName, 'power-off', args )
            elif self.LVtype == "RIGOL":
                self.RIGOLcmd(self.LVip, self.LVch, "turnOff()", args)
            elif self.LVtype == "CUSTOM":
                if args != None:
                    if args[0]=='cli':
                       self.conn.sendall( ( f'accepted: LV to set OFF.\n').encode() )
                self.logger.info( "CUSTOM LV command called: {}".format(self.LVOffCmd) )
                stdout = subprocess.run( self.LVOffCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT ).stdout.decode()
                self.logger.debug( stdout )
            else :
                self.conn.sendall( ('rejected: LV Power Off command  !\n').encode() ) 
    #--------------------------------------------------
    def setHV_On( self, args=None ):
        self.logger.debug( "Called setHV_On()" )
        self.checkCoverStatus()
        
        if self.isLV_On():
            if self.HVtype == "MPOD" :
                self.mpodcmd( self.mpodHVName, 'power-on', args )
            elif self.HVtype == "CUSTOM":
                if args != None:
                    if args[0]=='cli':
                       self.conn.sendall( ( f'accepted: HV to set ON.\n').encode() )
                self.logger.info( "CUSTOM HV command called: {}".format(self.HVOnCmd) )
                stdout = subprocess.run( self.HVOnCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT ).stdout.decode()
                self.logger.debug( stdout )
        elif args[0]=='cli':
            self.conn.sendall( ('rejected: LV is still OFF!\n').encode() )
            
        
        
    #--------------------------------------------------
    def setHV_Off( self, args=None ):
        self.logger.debug( "Called setHV_Off()" )
        if self.HVtype == "MPOD" :
            self.mpodcmd( self.mpodHVName, 'power-off', args )
        elif self.HVtype == "CUSTOM":
            if args != None:
               if args[0]=='cli':
                  self.conn.sendall( ( f'accepted: HV to set OFF.\n').encode() )
            self.logger.info( "CUSTOM HV command called: {}".format(self.HVOffCmd) )
            stdout = subprocess.run( self.HVOffCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT ).stdout.decode()
            self.logger.debug( stdout )
        
    #--------------------------------------------------
    def isHV_On( self, args=None ):
        self.logger.debug( "Called isHV_On()" )
        if self.HVtype == "MPOD" :
            out = self.snmpcmd( 'get', 'outputStatus', self.mpodHVChannel, '', '' )
            is_on = not ( out.find( '00 01' ) > 0 )
        elif self.HVtype == "CUSTOM":
            # HVIsOnCmd should return the measured voltage value. Non-zero value is expected when it's ON.
            
            self.logger.info( "CUSTOM HV command called: {}".format(self.HVIsOnCmd) )
            stdout = subprocess.run( self.HVIsOnCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT ).stdout.decode()
            self.logger.debug( "+++"+stdout+"+++" )
            return False if stdout.strip() == "0" else True

        return is_on

    #--------------------------------------------------
    def isLV_On( self, args=None ):
        self.logger.debug( "Called isLV_On()" )
        if self.LVtype == "MPOD" :
            out = self.snmpcmd( 'get', 'outputStatus', self.mpodLVChannel, '', '' )
            is_on = not ( out.find( '00 ' ) > 0 )
        elif self.LVtype == "RIGOL" :
            try:
                vol = float( self.RIGOLcmd(self.LVip, self.LVch, "measureVoltage()", args) )
                curr = float( self.RIGOLcmd(self.LVip, self.LVch, "measureCurrent()", args) )
                if vol > 1.0 and curr > 1.0 :
                    is_on = True
                else :
                    is_on = False
            except Exception as e:
                self.logger.error(f'Exception in RIGOL LV state checking\n\n: {str(e)}')
                self.sendEmail( 'RIGOL EXCEPTION', f"```\n{str(e)}\n```", category = 'ERROR' )
                return True
        elif self.LVtype == "CUSTOM":
            # LVIsOnCmd should return the measured current value. Non-zero value is expected when it's ON.
            
            self.logger.info( "CUSTOM LV command called: {}".format(self.LVIsOnCmd) )
            stdout = subprocess.run( self.LVIsOnCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT ).stdout.decode()
            self.logger.debug( "+++"+stdout+"+++" )
            return False if stdout.strip() == "0" else True

        return is_on
    #--------------------------------------------------
    def resetMPOD( self, args=None ):
        for group in [0, 64]:
            self.snmpcmd( 'set', 'groupsSwitch', group, 'i', 10 )
        
        if args and args[0]=='cli':
            self.conn.sendall( ('accepted: MPOD groups switches have been cleared.\n').encode() )

    #--------------------------------------------------
    def turnOffLV_afterInterval( self ):
        time.sleep(360)
        self.setLV_Off()
        
    #--------------------------------------------------
    def setVerbosity( self, args ):
        try:
            self.logger.setVerbosity( args[0] )
            self.conn.sendall( ('accepted: changed verbosity to {}.\n'.format( args[0] ) ).encode() )
        except:
            self.conn.sendall( ('rejected: invalid arg {}.\n'.format( args[0] ) ).encode() )

        
    #--------------------------------------------------------------------------------
    def sendEmail( self, title = '', message = '', category = 'INFO' ):
        # self.sendMessage( title = title, message = message, category = category, via = 'Email' )
        #self.sendMessage( title = title, message = message, category = category, via = 'Skype' )
        pass

    #--------------------------------------------------------------------------------
    def sendMessage( self, title = '', message = '', category = 'INFO' , via = 'Skype'):
        
        if via not in ['Skype','Email']:
            self.logger.error( 'Message cannot be sent via ' + via )

        try:
            base = os.getenv('PIDPATH')
            
            with open( '{}/config/Notification.json'.format( base ) ) as f:
                j = json.load( f )
                
            # cateogry: 'INFO', 'FATAL', 'ERROR', "WARNING'
            
            tstamp = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
            unix_begin=int( time.time()-600 )*1000
            unix_end=int( time.time()+600 )*1000
            grafana = f'[Grafana snapshot](http://192.168.100.104:3000/d/{self.grafanaID}/orgId=1&from={unix_begin}&to={unix_end}&refresh=10s)'
            messages = [ '\n', tstamp, message, grafana ]
            subject = '[ {} ]  {}  |  {}'.format( category, title, self.thisName )

            gmail = 'atlas.jp.itk.pixel.production@gmail.com'

            if via == 'Skype':
                msg = ' -- '.join( [ SkypeMsg.bold(subject) ] + messages )

                password = 'Thus-lcu7-5prA-8og&'
                hostname = os.uname()[1]
                chatid = j['chatid'][hostname]

                sk = Skype(gmail, password) # Tateyama ATLAS Japan ITk
                sk.chats[chatid].sendMsg(msg, rich=True)

            if via == 'Email':
                msg = MIMEText( '\n\n'.join( messages ) )

                password = 'nincgefvmkjgrcrs'
                to_email = ','.join( j['email_to'] )

                msg["Subject"] = subject
                msg["To"] = to_email
                msg["From"] = gmail

                server = smtplib.SMTP( 'smtp.gmail.com', 587 )
                server.set_debuglevel(False)
                server.ehlo()
                server.ehlo()
                server.starttls()
                server.login(gmail, password)
                server.sendmail( gmail, to_email, msg.as_string() )
                server.quit()
            
        except Exception as e:
            self.logger.error( pprint.pformat(traceback.format_exc()) )
            self.logger.error( 'Error in sending notification' )


    #--------------------------------------------------
    def sync(self, args):
        try:
            if args[0] == 'all':
                j = self.mon.toJson( args[0] )
            else:
                j = self.mon.toJson( float( args[0] ) )

            self.logger.debug( str(j) )
            self.conn.sendall( (str(j)+ '\n').encode()  )
        except:
            self.conn.sendall( ('sync error\n').encode() )
            self.logger.error( 'sync error' )
            
    
    #--------------------------------------------------
    def __init__(self, config):
        self.commands = { 'kill'             : self.kill,
                          'params'           : self.params,
                          'setT'             : self.setT,
                          'resetV'           : self.resetV,
                          'setKp'            : self.setKp,
                          'setKi'            : self.setKi,
                          'setKd'            : self.setKd,
                          'setMaxP'          : self.setMaxP,
                          'setMaxI'          : self.setMaxI,
                          'setMaxD'          : self.setMaxD,
                          'setIduration'     : self.setIduration,
                          'setTimeOffsetP'   : self.setTimeOffsetP,
                          'setTimeOffsetI'   : self.setTimeOffsetI,
                          'setTimeOffsetD'   : self.setTimeOffsetD,
                          'stopPID'          : self.stopPID,
                          'restartPID'       : self.restartPID,
                          'restoreInterlock' : self.restoreInterlock,
                          'unlock'           : self.unlock,
                          'setRelayOn'       : self.setRelayOn,
                          'setRelayOff'      : self.setRelayOff,
                          'resetDCS'         : self.resetDCS,
                          'setLV_On'         : self.setLV_On,
                          'setLV_Off'        : self.setLV_Off,
                          'setHV_On'         : self.setHV_On,
                          'setHV_Off'        : self.setHV_Off,
                          'resetMPOD'        : self.resetMPOD,
                          'setVerbosity'     : self.setVerbosity,
                          'sync'             : self.sync   }
        self.fatal = False
        self.kill = False
        self.jumpTolerance = 3
        self.V = 0.0
        self.Vprev = None
        self.setV = 0.0
        self.measV = 0.0
        self.logOutFile = 'tmp.log'
        self.isPidActivated = True
        self.isInterlockEnabled = False
        self.isInterlocked = False
        self.isInterlockedHighTemp = False
        self.isInterlockedHighDP   = False
        self.isInterlockedPeltierTripped   = False
        self.isInterlockedIrregularTargetTemp   = False
        self.lockStatus = -1
        self.is_dcsCtrl_good = True
        
        self.HVtype = config['HVtype']
        self.LVtype = config['LVtype']

        if self.HVtype == "MPOD" or self.LVtype == "MPOD":
            self.mpodIP = config['mpodIP']
            if self.HVtype == "MPOD":
                self.mpodHVChannel = config['mpodHV']
                self.mpodHVName = config['mpodHVname']
            if self.LVtype == "MPOD":
                self.mpodLVChannel = config['mpodLV']
                self.mpodLVName = config['mpodLVname']
        # self.HVip = config['HVip'] # Is it what?
        # self.HVch = config['HVch'] # Is it what?
        if self.LVtype == "RIGOL":
            self.LVip = config['LVip']
            self.LVch = config['LVch']
        if self.LVtype == "CUSTOM":
            self.LVOnCmd = config['LVOnCmd']
            self.LVOffCmd = config['LVOffCmd']
            self.LVIsOnCmd = config['LVIsOnCmd']
        if self.HVtype == "CUSTOM":
            self.HVOnCmd = config['HVOnCmd']
            self.HVOffCmd = config['HVOffCmd']
            self.HVIsOnCmd = config['HVIsOnCmd']

        self.thisName = config['name']
        self.grafanaID = config['grafanaID']


        self.peltError = False

        self.base             = os.getenv('PIDPATH')
        self.cvlccmd          = 'cvlc --play-and-exit '

        self.objName          = str( config['object'] )
        self.altName          = str( config['object_alt'] )
        self.minT             = float( config['minT'] )
        self.maxT             = float( config['maxT'] )
        self.T_default        = float( config['T_default'] )
        self.T_set            = self.T_default
        self.T_set_regulated  = self.T_default
        self.Kp               = float( config['Pcontrol'] )
        self.Ki               = float( config['Icontrol'] )
        self.Kd               = float( config['Dcontrol'] )
        self.interval         = float( config['interval'] )
        self.integralDuration = int( config['integralDuration'] )
        self.maxV             = float( config['maxV'] )
        self.maxPfeedback     = float( config['maxPfeedback'] )
        self.maxIfeedback     = float( config['maxIfeedback'] )
        self.maxDfeedback     = float( config['maxDfeedback'] )
        self.tOffsetP         = float( config['tOffsetP'] )
        self.tOffsetI         = float( config['tOffsetI'] )
        self.tOffsetD         = float( config['tOffsetD'] )
        self.interlockMargin  = float( config['interlockMargin'] )
        self.email_to         = [ 'ueyama.y.c90a@m.isct.ac.jp' ]
        

    #--------------------------------------------------
    def injectDependencies( self, mon, devices, logger, config ):
        self.mon         = mon
        self.peltCtrl    = devices['peltierControl']
        self.dcsCtrl     = devices['dcsControl']
        self.roomTempMon = devices['roomTemp']
        self.roomDPMon   = devices['roomDP']
        self.logger      = logger
        self.config      = config
    
        self.stabilityGaugingLength  = self.config['stability']['stabilityGaugingLength']
        self.deviationTolerance      = self.config['stability']['deviationTolerance']
        self.Tstdev                  = self.config['stability']['Tstdev']
        self.Vstdev                  = self.config['stability']['Vstdev']

    #--------------------------------------------------
    def calculateV(self, dV):
        
        self.mon.status = 'normal'

        self.Vprev = self.V
        self.V += dV
        
        if abs(self.V) > self.maxV:
            self.V = self.maxV * (-1.0 if self.V<0 else 1.0 )
            self.mon.status = 'overflow'
            
        try:
            if abs( self.mon.data['I'][-1] ) > 9.0:
                self.V = self.Vprev
        except:
            pass


    
    #--------------------------------------------------
    def calculateD( self, dTs, slope ):
        dT = abs( dTs[-1] )
        coeff = self.Kd

        if dT < 3.0:
            coeff *= max(1.0, 2.0 * dT - 1.0 )
        elif( dT > 3.0 ):
            coeff *= max(1.0, 15.0 / dT )
        
        return coeff * slope


    
    #--------------------------------------------------------------------------------
    def getExtrapolated( self, dTs, ts, offsetP, offsetI, offsetD ):

        span = min(max(60, 120 - int(abs(dTs[-1]))), len(ts))
        
        x = np.array( [ t - ts[-span] for t in ts[-span:] ] )
        y = np.array( KalmanSmoother( dTs[-span:], 0.05, 0.3 ).process() )
        weights = np.array( [ 1./0.25 for dT in dTs[-span:] ] )
        
        params = None
        
        if len(dTs) < 5:
            return 0, 0, 0
        else:
            params = np.polyfit( x, y, 2, w = weights )
            
        f = np.poly1d( params )
        d = np.polyder( f )

        i = reduce( add, [ f( x[-1] + dx ) for dx in list( range( int(offsetI-self.integralDuration), int(offsetI) ) ) ] ) / float( self.integralDuration )
        
        return f( x[-1] + offsetP ), i, d( x[-1] + offsetD )
    
    
    #--------------------------------------------------
    def shutdownMPOD( self, channel=None ):
        # temporary channels are hard-coded, to be fixed
        if channel == None:
            self.setHV_Off()
            self.setLV_Off()
            
        elif channel == 'HV':
            self.setHV_Off()
            
        elif channel == 'LV':
            self.setLV_Off('interlock')

    def increment( self, name ):


        # Exception: module NTC reading is disconnected.

        if self.isInterlockEnabled and self.mon.data[name][-1] is None:
            self.isInterlockedIrregularTargetTemp = True
            self.isInterlocked = True
            
            self.shutdownMPOD('HV')
            self.execDCS( partial(self.dcsCtrl.relayOff, 7 ) ) # HV
            self.execDCS( partial(self.dcsCtrl.relayOff, 6 ) ) # LV
            
            if not self.fatal:
                self.sendEmail( 'INTERLOCKED', 'ERROR: SW Interlock became active (Irregular Target Temp)!', category = 'ERROR' )
                os.environ["DISPLAY"]=":1"
                os.system( "${EQC_OPERATOR_DIR}/lib/Notifier.py --message \"ERROR::INTERLOCKED SW Interlock became active (Irregular Target Temp)!\" --color orange &" )
                os.system( "${EQC_OPERATOR_DIR}/gui/soundClient.py --host 192.168.100.104 --mode interlock &" )
                self.logger.error('SW Interlock became active (Irregular Target Temp)!')
            
            self.mon.data['isInterlockedIrregularTargetTemp'][-1] = int( self.isInterlockedIrregularTargetTemp )
            self.mon.data['isInterlockedPeltierTripped'][-1]      = int( self.isInterlockedPeltierTripped )
            self.mon.data['isInterlocked'                   ][-1] = int( self.isInterlocked )
            
            os.system( '{} {}/sounds/alarm1.mp3 &>/dev/null'.format( self.cvlccmd, self.base ) )
            P = 0.0
            I = 0.0
            D = 0.0

            dV = P + I + D
            self.logger.debug(f"w/ interlock {dV}={P}+{I}+{D}")
            
            self.calculateV( dV )
            
            for key, value in { 'T_set': float(self.T_set), 'T_set_regulated': float(self.T_set_regulated), 'Vp':P, 'Vi':I, 'Vd':D, 'dV' : dV, 'V': self.V }.items():
                if key in self.mon.data:
                    self.mon.data[key].append( value )
                else:
                    self.mon.data.update( { key: [value] } )
                    
            return

        
        #------------------------------------------
        
        ts = self.mon.data['elapsed'][-self.integralDuration:]
        Ts = [ self.T_set_regulated if abs( v ) > 50 else v for v in KalmanSmoother( self.mon.data[name][-self.integralDuration:], 0.02, 0.1 ).process() ]
        dTs = [ T - self.T_set_regulated for T in Ts ]
        
        #coeff = min( ( 0.2 + 0.8 * abs( self.T_set_regulated - self.mon.data['T_chiller'][-1] ) / 40.0 ), 1.0) if abs( self.mon.data[name][-1] - self.mon.data['T_chiller'][-1] ) < 10 else 1.0
        coeff = 1.0
        
        if False:
            val_extrap, integral_extrap, slope_extrap = self.getExtrapolated( dTs, ts, self.tOffsetP, self.tOffsetI, self.tOffsetD )

            P = coeff * self.Kp * dTs[-1]
            I = coeff * self.Ki * integral_extrap
            D = coeff * self.calculateD( dTs, slope_extrap )
            
        else:

            val_extrap, integral_extrap, slope_extrap = self.getExtrapolated( dTs, ts, self.tOffsetP, self.tOffsetI, self.tOffsetD )
            integral = reduce( add, dTs ) / float( len(dTs) ) * (1.0 if abs( dTs[-1] ) < 5.0 else 5.0/abs( dTs[-1] ) )
            slope    = getSlope( dTs, ts )
            
            P = coeff * self.Kp * val_extrap
            I = coeff * self.Ki * integral_extrap
            D = coeff * self.calculateD( dTs, slope )

        # regulate the I-feedback (diminish at large deltaT)
        I = I * math.exp( - dTs[-1]*dTs[-1]/4.0 )

        P = min( P, self.maxPfeedback ) if P > 0 else -1.0 * min( abs(P), self.maxPfeedback )
        I = min( I, self.maxIfeedback ) if I > 0 else -1.0 * min( abs(I), self.maxIfeedback )
        D = min( D, self.maxDfeedback ) if D > 0 else -1.0 * min( abs(D), self.maxDfeedback )

        if not self.isPidActivated:
            P = 0.0
            I = 0.0
            D = 0.0

        dV = P + I + D
            
        self.logger.debug(f"w/o interlock {dV}={P}+{I}+{D}")
        self.calculateV( dV )

        for key, value in { 'T_set': float(self.T_set), 'T_set_regulated': float(self.T_set_regulated), 'Vp':P, 'Vi':I, 'Vd':D, 'dV' : dV, 'V': self.V, 'slope':slope }.items():
            if key in self.mon.data:
                self.mon.data[key].append( value )
            else:
                self.mon.data.update( { key: [value] } )
                self.logger.debug(f'{key}: {[value]}')
                
        
    #--------------------------------------------------
    def isHumiditySafe(self):
        if self.mon.data[self.objName][-1] is not None and self.mon.data[self.objName][-1] > -99:
            return self.mon.data['DP_carrier'][-1] < self.mon.data[self.objName][-1] - self.interlockMargin
        else:
            return True
        
    #--------------------------------------------------
    def isTemperatureSafe(self):
        if self.isInterlockEnabled:
            try: 
                self.logger.debug(f"temp w/ interlock: {self.mon.data[self.objName][-1]}")
                return self.mon.data[self.objName][-1] < 40.0 and self.mon.data['T_chuck'][-1] < 40.0 and self.mon.data['T_sink'][-1] < 40.0
            except:
                return True
        else:
            try:
                self.logger.debug(f"temp w/o interlock: {self.mon.data[self.altName][-1]}, {self.mon.data['T_sink'][-1]}")
                return self.mon.data[self.altName][-1] < 40.0 and self.mon.data['T_sink'][-1] < 40.0
            except:
                return True

    
    #--------------------------------------------------
    def checkCoverStatus(self):
        
        status_prev = self.lockStatus
        
        self.lockStatus = self.execDCS( self.dcsCtrl.lockStatus )
        
        if self.lockStatus == 1: # Cover is closed
            self.isInterlockEnabled = True

            if self.lockStatus != status_prev:
                self.sendEmail( 'Cover has been closed', category = 'INFO' )
            
        else:
            if self.lockStatus != status_prev:
                self.sendEmail( 'Cover is open', category = 'WARNING' )
                
            self.logger.debug("pidLoop(): Cover is open --> no interlock")
            self.isInterlockEnabled = False
                
        
    #--------------------------------------------------
    def recoverPeltierControl(self):
        self.logger.info( 'Detected peltierControl device' )
        time.sleep(30)
        
        self.peltCtrl = TexioPFRControl( { "port": self.config['devices']['peltierControl']['port'],
                                           "address": self.config['devices']['peltierControl']['address'],
                                           "baud":self.config['devices']['peltierControl']['baud'] } )
        
        time.sleep(1)
        self.peltCtrl.resetAlarms()
        time.sleep(3)
        
        status = self.peltCtrl.getStatus()
        self.logger.verbose( 'pidLoop(): peltier controller status = ' + str( status ) )
        
        status = self.peltCtrl.setV(0.0)
        self.peltCtrl.setOff()
        time.sleep(1)
        self.peltCtrl.setOn()
        time.sleep(3)
        status = self.peltCtrl.setV( self.V )
        
        self.logger.info( 'Succeeded in re-invoking peltierControl' )
        self.sendEmail( 'Succeeded in re-invoking peltierControl' )
        
        self.peltError = False
        
        return
        
        
    #--------------------------------------------------
    def pidLoop(self):
        
        self.logger.info('pidLoop(): started')
        
        self.checkCoverStatus()
        
        self.isInterlocked = False
        self.isInterlockedHighTemp = False
        self.isInterlockedHighDP = False
        self.isInterlockedPeltierTripped = False
        self.isPidActivated = True


        entities = [ 'elapsed', 'T_module', 'T_chuck', 'T_case', 'T_carrier', 'T_chiller', 'T_sink', 
                     'RH_carrier', 'DP_carrier',
                     'T_chuck_ref', 'T_case_ref', 'T_chiller_ref', 'T_sink_ref', 
                     'interlockStatus', 'RS_heater', 'RS_peltier', 'RS_chiller', 'RS_pelPlus', 'RS_pelMinus',
                     'RS_lowVoltage', 'RS_highVoltage', 'RS_coolingBoxLock', 'RS_lockState1', 'RS_lockState2', 'RS_chillerAlert' ]
        
        while True:

            #--------------------------------------------------
            # Acquire new monitoring data
            #--------------------------------------------------
            
            self.logger.debug( " pidLoop(): Acquire new monitoring data..." )

            try:
                elapsed = self.mon.getElapsed()
                
                module_temp = -9999
                if self.dcsCtrl.m_temps['module'] is not None and self.dcsCtrl.m_temps['module'] > -60:
                    module_temp = self.dcsCtrl.m_temps['module']
                
                self.execDCS( self.dcsCtrl.readAll )
                updates = [ elapsed,
                            module_temp,
                            self.dcsCtrl.m_temps['head'],
                            self.dcsCtrl.m_temps['case'],
                            self.dcsCtrl.m_temps['carrier'],
                            self.dcsCtrl.m_temps['chiller'],
                            self.dcsCtrl.m_temps['sink'],
                            self.dcsCtrl.m_rh['carrier'],
                            self.dcsCtrl.m_dp['carrier'],
                            self.dcsCtrl.m_temps_ref['head'],
                            self.dcsCtrl.m_temps_ref['case'],
                            self.dcsCtrl.m_temps_ref['chiller'],
                            self.dcsCtrl.m_temps_ref['sink'],
                            self.dcsCtrl.m_interlockStatus,
                            self.dcsCtrl.m_rState['heater'],
                            self.dcsCtrl.m_rState['peltier'],
                            self.dcsCtrl.m_rState['chiller'],
                            self.dcsCtrl.m_rState['pelPlus'],
                            self.dcsCtrl.m_rState['pelMinus'],
                            self.dcsCtrl.m_rState['lowVoltage'],
                            self.dcsCtrl.m_rState['highVoltage'],
                            self.dcsCtrl.m_rState['coolingBoxLock'],
                            self.dcsCtrl.m_rState['lockState1'],
                            self.dcsCtrl.m_rState['lockState2'],
                            self.dcsCtrl.m_rState['chillerAlert']
                ]

            except:
                if not self.fatal:
                    self.sendEmail( 'Exception raised in increment()', traceback.format_exc(), category='FATAL'  )
                    self.logger.error('Exception raised in increment()')
                    self.logger.error( pprint.pformat( traceback.format_exc() ) )
                    
                self.fatal = True
                
                time.sleep(3)
                continue
                
            for e, u in zip(entities, updates):
                if e in self.mon.data.keys():
                    # temporary protection for the chiller temp jump
                    if e == 'T_chiller' and abs( self.mon.data[e][-1] - u  ) > 4.0 :
                        self.mon.data[e].append( self.mon.data[e][-1] )
                    else:
                        self.mon.data[e].append( u )
                    
                else:
                    self.mon.data.update( { e:[u] } )

            # regulate the actual target value of PID feedback
            
            self.logger.debug( "pidLoop(): regulate the actual target value of PID feedback" )
            if self.T_set < self.mon.data['DP_carrier'][-1]+self.interlockMargin+2:
                self.T_set_regulated = self.mon.data['DP_carrier'][-1]+self.interlockMargin+2
            elif self.T_set > self.maxT:
                self.T_set_regulated = self.maxT
            elif self.T_set < self.minT:
                self.T_set_regulated = self.minT
            else:
                self.T_set_regulated = self.T_set
                
            #--------------------------------------------------
            # Judge HW interlock
            #--------------------------------------------------
            
            if len( self.mon.data['interlockStatus'] ) > 1:
                status_prev = self.mon.data['interlockStatus'][-2]
                status_curr = self.mon.data['interlockStatus'][-1]
                if status_curr > 0 and status_prev < status_curr :
                    
                    status_map = ['INTERLOCK_NORMAL',
                                  'INTERLOCK_HIGHTEMP_LV',
                                  'INTERLOCK_HIGHTEMP_PELTIER',
                                  'INTERLOCK_HIGHTEMP_IDLE',
                                  'INTERLOCK_HIGHDEW_PELCHILLER',
                                  'INTERLOCK_HIGHDEW_IDLE',
                                  'INTERLOCK_CHILLER_ALERT',
                                  'INTERLOCK_NTC_DISONN' ]

                    self.logger.error("HW Interlock became active: mode = {} status_code = {}".format( status_curr, status_map[status_curr] ) )
                    self.isPidActivated = False
                    self.sendEmail( 'HARDWARE INTERLOCK', 'HW Interlock became active!\nMode = {} status_code = {}'.format( status_curr, status_map[status_curr] ), category = 'ERROR' )
                    os.system( "${EQC_OPERATOR_DIR}/gui/soundClient.py --host 192.168.100.104 --mode interlock &" )
                    os.environ["DISPLAY"]=":1"
                    os.system( f"${{EQC_OPERATOR_DIR}}/lib/Notifier.py --message \"ERROR::HARDWARE_INTERLOCK HW Interlock became active! Mode = {status_curr} status_code = {status_map[status_curr]}\" --color red &" )
                    
                    # Protect Module from heating up (redundant)
                    if status_map[status_curr] in [ 'INTERLOCK_CHILLER_ALERT', 'INTERLOCK_HIGHTEMP_LV', 'INTERLOCK_HIGHTEMP_PELTIER', 'INTERLOCK_HIGHTEMP_IDLE' ]:
                        self.shutdownMPOD()
                        self.execDCS( partial(self.dcsCtrl.relayOff, 7 ) ) # HV
                        self.execDCS( partial(self.dcsCtrl.relayOff, 6 ) ) # LV
                    if status_map[status_curr] in [ 'INTERLOCK_HIGHDEW_PELCHILLER', 'INTERLOCK_HIGHDEW_IDLE' ]:
                        self.shutdownMPOD('HV')
                        self.execDCS( partial(self.dcsCtrl.relayOff, 7 ) ) # HV
                        time.sleep(300)
                        self.shutdownMPOD('LV')
                        self.execDCS( partial(self.dcsCtrl.relayOff, 6 ) ) # LV
            

            #--------------------------------------------------
            # Judge SW interlock
            #--------------------------------------------------

            # check the cover switch status
            # once the cover is closed, interlock must be enabled.

            self.checkCoverStatus()
            
            if self.isInterlockEnabled:
                if self.isInterlocked:
                    pass
                else:
                    self.isInterlockedHighTemp = not self.isTemperatureSafe()
                    self.isInterlockedHighDP   = not self.isHumiditySafe()
                    
                    try:
                        self.isInterlockedPeltierTripped = self.peltCtrl.isTripped()
                    except:
                        self.isInterlockedPeltierTripped = True
                    
                    
                    self.isInterlocked = self.isInterlockedHighTemp or self.isInterlockedHighDP or self.isInterlockedPeltierTripped or not self.is_dcsCtrl_good
                    
                    # Force shutdown LV
                    if self.isInterlockedHighTemp:
                        
                        tHV = threading.Thread( target = self.setHV_Off, args = () )
                        tHV.start()
                        tLV = threading.Thread( target = self.setLV_Off, args = ('interlock') )
                        tLV.start()
                        self.isPidActivated = False
                        self.V = 0.0
                            
                        self.logger.error("LV and HV turned OFF!")
                        self.logger.warning("SW Interlock became active!")
                        self.sendEmail( 'INTERLOCKED', 'ERROR: SW Interlock became active! (High Temp)\n\nLV and HV turned OFF', category = 'ERROR' )
                        os.environ["DISPLAY"]=":1"
                        os.system( f"${{EQC_OPERATOR_DIR}}/lib/Notifier.py --message \"ERROR::INTERLOCKED SW Interlock became active! (High Temp) LV and HV turned OFF\" --color orange &" )
                        os.system( "${EQC_OPERATOR_DIR}/gui/soundClient.py --host 192.168.100.104 --mode interlock &" )
                        
                    if self.isInterlockedHighDP:
                        
                        tHV = threading.Thread( target = self.setHV_Off, args = () )
                        tHV.start()
                        
                        self.T_set = self.T_set + 5.0
                        
                        if self.mon.data[self.objName][-1] > 20. :
                            thread = threading.Thread( target = self.turnOffLV_afterInterval, args = () )
                            thread.start()
                        
                        self.logger.warning("SW Interlock became active!")
                        self.sendEmail( 'INTERLOCKED', 'ERROR: SW Interlock became active! (High DP)\n\nHV turned OFF', category = 'ERROR' )
                        os.environ["DISPLAY"]=":1"
                        os.system( f"${{EQC_OPERATOR_DIR}}/lib/Notifier.py --message \"ERROR::INTERLOCKED SW Interlock became active! (High DP) HV turned OFF\" --color orange &" )
                        os.system( "${EQC_OPERATOR_DIR}/gui/soundClient.py --host 192.168.100.104 --mode interlock &" )
                    
                    if self.isInterlockedPeltierTripped:
                        tHV = threading.Thread( target = self.setHV_Off, args = () )
                        tHV.start()
                        
                        tLV = threading.Thread( target = self.setLV_Off, args = ('interlock') )
                        tLV.start()
                        
                        self.isPidActivated = False
                        self.V = 0.0
                        
                        self.logger.error("Shut down LV Relay!")
                        self.logger.warning("SW Interlock became active!")
                        self.sendEmail( 'INTERLOCKED', 'ERROR: SW Interlock became active! (Peltier Tripped)\n\nLV and HV turned OFF', category = 'ERROR' )
                        os.environ["DISPLAY"]=":1"
                        os.system( f"${{EQC_OPERATOR_DIR}}/lib/Notifier.py --message \"ERROR::INTERLOCKED SW Interlock became active! (Peltier Tripped) LV and HV turned OFF\" --color orange &" )
                        os.system( "${EQC_OPERATOR_DIR}/gui/soundClient.py --host 192.168.100.104 --mode interlock &" )
                        
            else:
                self.logger.debug( "pidLoop(): interlock inactive --> force setT at {} degC.".format( self.T_default ) )
                
                
                self.isInterlockedHighDP   = False # Not applied when the cover is open
                self.isInterlockedHighTemp = not self.isTemperatureSafe()
                
                if self.isInterlockedHighTemp:
                    tHV = threading.Thread( target = self.setHV_Off, args = () )
                    tHV.start()
                    
                    tLV = threading.Thread( target = self.setLV_Off, args = ('interlock') )
                    tLV.start()
                    
                    self.V = 0.0
                
                    self.logger.error("Shut down LV and HV Relays! Peltier voltage is set zero!")
                    self.logger.warning("HW Interlock [HighTemp] became active!")
                    
                self.isInterlocked = self.isInterlockedHighTemp
                self.T_set = self.T_default
                self.T_set_regulated = self.T_default


            # add interlock and PID activation status to mon.data
            for var in [ 'coverStatus', 'isInterlockEnabled', 'isInterlocked', 'isInterlockedHighTemp', 'isInterlockedHighDP', 'isInterlockedIrregularTargetTemp', 'isInterlockedPeltierTripped', 'isPidActivated' ]:
                if not var in self.mon.data:
                    self.mon.data.update( { var: [] } )
            
            self.mon.data['coverStatus']                      .append( self.lockStatus )
            self.mon.data['isInterlockEnabled']               .append( int(self.isInterlockEnabled) )
            self.mon.data['isInterlocked']                    .append( int(self.isInterlocked) )
            self.mon.data['isInterlockedHighTemp']            .append( int(self.isInterlockedHighTemp) )
            self.mon.data['isInterlockedHighDP']              .append( int(self.isInterlockedHighDP) )
            self.mon.data['isInterlockedIrregularTargetTemp'] .append( int(self.isInterlockedIrregularTargetTemp) )
            self.mon.data['isInterlockedPeltierTripped']      .append( int(self.isInterlockedPeltierTripped) )
            self.mon.data['isPidActivated']                   .append( int(self.isPidActivated) )
            
            self.mon.data['RS_peltier'][-1] = int( self.mon.data['RS_chillerAlert'][-1] or not self.peltError )
            
            # First check the status of the power supply
            try:
                
                if not self.peltError:
                    status = self.peltCtrl.getStatus()
                    
            except:
                self.logger.error( pprint.pformat( e ) )
                self.peltCtrl.m_ctrl.close()
                
                if not self.peltError:
                    self.sendEmail( 'Exception raised in increment()', traceback.format_exc(), category='FATAL'  )
                    self.logger.error('Exception raised in pidLoop() after peltCtrl.getStatus()')
                    self.peltError = True
                    self.peltCtrl.m_ctrl.close()
                self.shutdownMPOD()
                
            if self.peltError:
                
                try:
                    
                    if os.path.exists( self.config['devices']['peltierControl']['port'] ):
                        self.recoverPeltierControl()
            
                except Exception as e:
                    self.logger.error( pprint.pformat( e ) )
                    self.logger.error( pprint.pformat( traceback.format_exc() ) )
                    self.logger.error( 'Exception raised in re-invoking peltierControl!' )
                    self.isInterlocked = True
                    self.shutdownMPOD()

                    self.peltError = True
                        
                    database    = self.config['influxDB']['database']
                    hostname    = self.config['influxDB']['hostname']
                    port        = self.config['influxDB']['port']
                    measurement = self.config['influxDB']['measurement']
                    
                    self.mon.recordToDB( self.logger,
                                         host = hostname,
                                         port = port,
                                         database = database,
                                         measurement = measurement )
                    
                    continue

            #--------------------------------------------------
            # Case-1: interlock is active
            #--------------------------------------------------
            if self.isInterlocked:

                self.logger.debug( "pidLoop(): interlock inactive --> regulate the setT_actual" )
                self.T_set_regulated = min( max( self.mon.data['DP_carrier'][-1]+self.interlockMargin+5, self.T_set ), self.maxT )
                
                try:
                    self.increment( self.objName )
                except:
                    self.logger.error( 'unexpected error raised. you may want to check /tmp/pidServer_console.txt' )
                    print(traceback.format_exc())
                    if not self.fatal:
                        self.sendEmail( 'Exception raised in increment()', traceback.format_exc(), category='FATAL'  )
                        self.logger.error('Exception raised in pidLoop() after increment()')
                    self.fatal = True
                    time.sleep(3)
                    continue
                
            #--------------------------------------------------
            # Case-2: interlock is inactive, and PID is active (normal mode)
            #--------------------------------------------------
            else:
                
                try:
                    # if interlock is active, PID looks at the target object temperature
                    # otherwise, PID looks at the vacuum chuck temperature (i.e. the case is open)
                    if self.isInterlockEnabled:
                        self.increment( self.objName )
                    else:
                        self.logger.debug( "pidLoop(): interlock inactive: incretement PID using vac.chuck" )
                        self.increment( self.altName )
                        self.logger.debug("test13")
                except:
                    if not self.fatal:
                        self.sendEmail( 'Exception raised in increment()', traceback.format_exc(), category='FATAL'  )
                        self.logger.error('Exception raised in pidLoop() after increment()')
                    self.logger.error( 'unexpected error raised. you may want to check /tmp/pidServer_console.txt' )
                    print(traceback.format_exc())
                    self.fatal = True
                    time.sleep(3)
                    continue

            
            #--------------------------------------------------
            # Execute feedback to the power supply and relay switches
            #--------------------------------------------------
            
            self.logger.debug( "pidLoop(): reflect feedback result to peltier control" )
            # Flip the polarity of the peltier current
            if( self.Vprev < 0.0 and self.V >= 0 ):
                self.logger.debug("test11")
                if not self.isInterlockedPeltierTripped:
                    self.peltCtrl.setV( 0.0 )
                    time.sleep(1)
                    self.execDCS( self.dcsCtrl.pelPolNormal )
                    self.logger.debug("test12")
            elif( self.Vprev >= 0.0 and self.V < 0 ):
                if not self.isInterlockedPeltierTripped:
                    self.peltCtrl.setV( 0.0 )
                    time.sleep(1)
                    self.execDCS( self.dcsCtrl.pelPolInverted )
                    
            # Cross-check the polarity of the peltier voltage
            self.execDCS( self.dcsCtrl.relayStates )
            if self.V >= 0:
                if not ( self.dcsCtrl.m_rState['pelPlus'] == 0 and self.dcsCtrl.m_rState['pelMinus'] == 0 ):
                    self.peltCtrl.setV( 0.0 )
                    time.sleep(1)
                    self.execDCS( self.dcsCtrl.pelPolNormal )
            else:
                if not ( self.dcsCtrl.m_rState['pelPlus'] == 1 and self.dcsCtrl.m_rState['pelMinus'] == 1 ):
                    self.peltCtrl.setV( 0.0 )
                    time.sleep(1)
                    self.execDCS( self.dcsCtrl.pelPolInverted )

                    
            time.sleep(0.2)
            
            # Finally set the next voltage value
            try:
                if not self.isInterlockedPeltierTripped:
                    self.logger.debug(f"test14:{self.V}")
                    self.peltCtrl.setV( abs( self.V ) )
                    self.logger.debug(f"test15:{self.V}")
            except Exception as e:
                self.logger.error( pprint.pformat( e ) )
                self.logger.error( 'Failed in peltCtrl.setV()' )
            
            time.sleep(0.2)
            
            # Measure the current and update the data
            self.logger.debug( "pidLoop(): update monitoring" )
            I = 0.0
            try:
                if not self.isInterlockedPeltierTripped:
                   # I = self.peltCtrl.getI()
                    I = self.peltCtrl.getI()['meas']
                    if( self.V < 0 ): I = -I
                    
                    Vinfo = self.peltCtrl.getV()
                    self.setV = Vinfo.get('set')
                    self.measV = Vinfo.get('meas')
                    if( self.V < 0 ):
                        self.setV = -self.setV
                        self.measV = -self.measV
                    
            except Exception as e:
                self.logger.error( pprint.pformat( e ) )
                self.logger.error( 'Failed in peltCtrl.getI()' )
                
                
            if 'I' in self.mon.data:
                self.mon.data['I'].append( I )
            else:
                try:
                    self.mon.data.update( {'I':[ I ] } )
                except:
                    self.mon.data.update( {'I':[ None ] } )

            if 'setV' in self.mon.data:
                self.mon.data['setV'].append( self.setV )
            else:
                try:
                    self.mon.data.update( {'setV':[ self.setV ] } )
                except:
                    self.mon.data.update( {'setV':[ None ] } )

            if 'measV' in self.mon.data:
                self.mon.data['measV'].append( self.measV )
            else:
                try:
                    self.mon.data.update( {'measV':[ self.measV ] } )
                except:
                    self.mon.data.update( {'measV':[ None ] } )
                    
            if abs( self.V - self.measV ) > 2.0:
                self.isInterlockedPeltierTripped = True

            # Evaluate the stability
            self.logger.debug( "pidLoop(): stability calculation" )
            
            trend = None
            if self.isInterlockEnabled:
                trend      = np.array( self.mon.data[self.objName][-self.stabilityGaugingLength:], dtype=float )
            else:
                trend      = np.array( self.mon.data[self.altName][-self.stabilityGaugingLength:], dtype=float )
                
            trendMean  = trend.mean()
            TtrendStdev = trend.std()
            VtrendStdev = np.array( self.mon.data['V'][-self.stabilityGaugingLength:] ).std()
            
            try:
                slopeTrend = np.array( self.mon.data['slope'][-self.stabilityGaugingLength:] ).mean()
            except Exception as e:
                self.logger.warning( str(e) )
                slopeTrend = 0.0
        
            self.isConverging   = 1 if ( len(trend) >= self.stabilityGaugingLength and
                                         abs(trendMean - self.T_set) < self.deviationTolerance and
                                         TtrendStdev < self.Tstdev and abs(slopeTrend) < 0.025 ) else 0

            if self.T_set_regulated != self.T_set:
                self.isConverging = False

            if 'isConverging' in self.mon.data:
                self.mon.data['T_rms'].append( TtrendStdev )
                self.mon.data['V_rms'].append( TtrendStdev )
                self.mon.data['isConverging'].append( int( self.isConverging ) )

                if len( self.mon.data['isConverging'] ) >= 80:
                    tail = self.mon.data['isConverging'][-80:]
                    convergingsCount = len( [ k for k in tail if k == 1 ] )
                    isStable = 1 if ( convergingsCount > 60 ) else 0
                else:
                    isStable = 0
                    
                self.mon.data['isStable'].append( isStable )

                    
            else:
                self.mon.data.update( {'T_rms': [ TtrendStdev ] } )
                self.mon.data.update( {'V_rms': [ VtrendStdev ] } )
                self.mon.data.update( {'isConverging':[ 0 ] } )
                self.mon.data.update( {'isStable': [ 0 ] } )
            
            
            #--------------------------------------------------
            # Heater control
            #--------------------------------------------------

            if (self.T_set_regulated > 28) and (not self.isInterlocked) and (self.T_set_regulated > self.mon.data['T_case'][-1] + 5):
                
                self.logger.debug( f'T_set_regulated = {self.T_set_regulated}, T_case = {self.mon.data["T_case"][-1]} --> heater on' )
                self.execDCS( partial( self.dcsCtrl.relayOn, 1 ) ) # Heater
                
            else:
                
                self.execDCS( partial( self.dcsCtrl.relayOff, 1 ) ) # Heater


            # Write records to InfluxDB
            self.logger.debug( "pidLoop(): logging" )
            self.mon.recordToDB( self.logger,
                                 host = self.config['influxDB']['hostname'],
                                 port = self.config['influxDB']['port'],
                                 database = self.config['influxDB']['database'],
                                 measurement = self.config['influxDB']['measurement'] )
            
            # Flow control with signals
            if self.kill:
                self.sendEmail( 'SERVER_KILLED', 'WARNING: PidServer Killed!', category = 'WARNING' )
                self.logger.warning( 'pidLoop(): killed' )
                return
            
            if self.fatal:
                self.logger.error( 'pidLoop(): fatal error' )
            
            self.logger.debug(f'{self.config['influxDB']['measurement']}') 
            # Wait for next round
            wait = self.interval - ( self.mon.getElapsed() - self.mon.data['elapsed'][-1] )
            if( wait < 0 ) : wait = self.interval
            time.sleep( wait )


    #--------------------------------------------------
    def userCommand(self):
        
        while True:
            tokens = str(  raw_input() ).split(' ')
            cmdName = tokens[0]
            
            if cmdName in self.commands.keys():
                command = self.commands[ cmdName ]
                args    = tokens[1:]
                command( args )
            else:
                self.logger.warning( 'invalid control command' )
                continue
            
            if self.kill == True:
                break
            
            if self.fatal == True:
                break
    
    #--------------------------------------------------
    def listenCommands(self, port = 50007):

        while True:
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.bind( ('', port) )
                s.listen(1)
            except:
                time.sleep(3)
                continue
            
            
            while True:
                self.conn = ( s.accept() )[0]
                self.logger.debug( 'client has connected' )
                while True:

                    try:
                        message = self.conn.recv(1024).decode()
                        tokens  = message.split(' ')
                        cmdName = tokens[0]
                        
                        if cmdName in self.commands.keys():
                            command = self.commands[ cmdName ]
                            args    = tokens[1:]

                            if not( cmdName == 'sync' or cmdName == 'params' ):
                                self.logger.info( 'received command from client: ' + message )
                                
                            command( args )
                                
                        elif cmdName.find(".q") == 0:
                            self.conn.sendall( 'accepted client quit\n'.encode() )
                            self.conn.close()
                            self.logger.debug( 'client has disconnected' )
                            
                        else:
                            self.conn.sendall( ('invalid control command: ' + message + '\n').encode() )
                            continue
                        
                        if self.kill == True:
                            time.sleep(3)
                            self.logger.warning( 'kill was signaled from the client' )
                            self.conn.close()
                            self.logger.debug( 'client has disconnected' )
                            s.close()
                            self.logger.debug( 'listening closed' )
                            return

                    except socket.error:
                        self.logger.debug( 'socket error --> recovery' )
                        break

                    except:
                        print(traceback.format_exc()) 
                        self.logger.error( 'unexpected error raised. you may want to check /tmp/pidServer_console.txt' )
                        self.sendEmail( 'unexpected error raised. you may want to check /tmp/pidServer_console.txt', traceback.format_exc(), category = 'FATAL' )
                        break
