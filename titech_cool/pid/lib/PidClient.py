import sys
import optparse
import traceback
from cmd import Cmd

import socket
import threading
import time
import re
import json
import psutil

from PidControl import *

EQC_OPERATOR_DIR = os.environ.get('EQC_OPERATOR_DIR', '')
INSTITUTION = os.environ.get('INSTITUTION', '')
HOSTNAME = os.environ.get('HOSTNAME', '')



class PidClient(Cmd):
    intro  = 'Press Ctrl-D to close the prompt.\n Type help to glance guides.'
    prompt = 'CoolingBox >>> '
    def __init__(self, logger, mon, options, debug = True):
        Cmd.__init__(self, completekey='tab')
        self.kill = False
        self.debug = debug
        self.logger = logger
        self.mon    = mon
        self.isBusy = False
        self.params = {}
        self.options = options

    def socketCall( self, message, writeOut = True ):

        while True:
            if self.isBusy:
                print( 'waiting for other busy thread...' )
                time.sleep(0.5)
            else:
                break

        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect( (self.options.host, self.options.port) )
        
        if not 'params' in message:
            self.logger.info( 'sending message to the server: ' + message )
        
        self.isBusy = True
        sock.sendall( message.encode() )
        
        buf = b''
        while True:
            data = sock.recv(1024)
            self.logger.debug( 'received packet from server: ' + repr(data) )
            buf += data
            if repr(data).find( '\\n' ) < 0:
                continue
            else:
                break
            
        msg = buf.decode("utf-8").replace('\\n', '').replace("'", "")
        
        if writeOut:
            if not 'params' in message:
                print( '[server]: ' + msg )
            else:
                j = json.loads( msg )
                print( json.dumps( j, sort_keys = False, indent = 4 ) )
                print( '{' )
                for key, d in sorted( self.mon.data.items() ):
                    print( '    "' + key + '": ' + str( d[-1] ) )
                print( '}' )
        
        sock.close()
        
        self.isBusy = False
        
        return msg
        

    def debugMsg( self, name, arg ):
        pass #if self.debug: print( name + "() called; arg = " + arg )

    def do_kill_server(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        choice = input( '\n******************\nAre you really sure to kill the server?\nPlease remember you need to assure the safety of the module after killing the PID service!!!\n******************\n\nChoose [y/N]: ').lower()
        if choice.find('y') == 0:
            self.kill = True
            self.socketCall( 'kill' )
            return True
        else:
            print('\nCancelled kill_server\n')

    def do_relaunch_server(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        choice = input( '\n******************\nAre you really sure to kill and the ongoing server and relaunch a new one?\n******************\n\nChoose [y/N]: ').lower()
        if choice.find('y') == 0:
            
            self.socketCall( 'kill' )
            
            self.isBusy = True
            
            print('... Please wait for server to completely shutdown...')
            
            while True:
                
                isRunning = False
                for proc in psutil.process_iter():
                    if 'pidServer' in proc.name():
                        isRunning = True
                        break
                
                if isRunning:
                    time.sleep(1)
                    continue
                
                break
                
            time.sleep(5)
            print('... Launching a new server...')
            os.system( 'pidServer.py ' + arg )
            
            print('... Please wait for server to completely launch...')
            time.sleep(15)
            self.isBusy = False
            
            return True
            
        else:
            print('\nCancelled relaunch_server\n')

    def do_quit(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        self.socketCall( '.q' )
        self.kill = True
        return True
        
    def do_EOF(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        print ('')
        if not self.kill:
            self.socketCall( '.q' )
            self.kill = True
        return True

    def do_exit(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        print ('')
        if not self.kill:
            self.socketCall( '.q' )
            self.kill = True
        return True

    def do_log(self, arg):
        
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.mon.printOut( int( args[0] ) )
        except:
            self.mon.printOut()

    def do_params(self, arg = ''):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            msg = self.socketCall( 'params', writeOut = False if arg == 'quiet' else True )
            self.params = json.loads( msg )
        except:
            print(traceback.format_exc())
            self.help_params()

    def do_setT(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            if args[0] == 'cold':
                args[0] = '-15'
            elif args[0] == 'warm':
                args[0] = '20'
            elif args[0] == 'exchange':
                args[0] = '30'
            elif args[0] == 'cold_startup':
                args[0] = '-35'
                
            self.socketCall( 'setT ' + args[0] )
        except:
            self.help_setT()

    def do_setKp(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setKp ' + args[0] )
        except:
            self.help_setKp()

    def do_setKi(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setKi ' + args[0] )
        except:
            self.help_setKi()

    def do_setKd(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setKd ' + args[0] )
        except:
            self.help_setKd()

    def do_setMaxP(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setMaxP ' + args[0] )
        except:
            self.help_setMaxP()

    def do_setMaxI(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setMaxI ' + args[0] )
        except:
            self.help_setMaxI()

    def do_setMaxD(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setMaxD ' + args[0] )
        except:
            self.help_setMaxD()

    def do_setIduration(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setIduration ' + args[0] )
        except:
            self.help_setIduration()

    def do_setTimeOffsetP(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setTimeOffsetP ' + args[0] )
        except:
            self.help_setTimeOffsetP()

    def do_setTimeOffsetI(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setTimeOffsetI ' + args[0] )
        except:
            self.help_setTimeOffsetI()

    def do_setTimeOffsetD(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'setTimeOffsetD ' + args[0] )
        except:
            self.help_setTimeOffsetD()

    def do_resetV(self, arg):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            args = arg.split()
            self.socketCall( 'resetV ' + args[0] )
        except:
            self.help_setT()

    def do_unlock( self, arg ):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            msg = self.socketCall( 'unlock' )
            
            if msg.find('accepted') >= 0:
                subprocess.run(f'{EQC_OPERATOR_DIR}/lib/serialNumberSwitcher.py', shell=True)
                
            return msg
        except:
            self.help_unlock()
            
        
    def do_resetDCS(self, arg):
        try:
            args = arg.split()
            self.socketCall( 'resetDCS' )
        except:
            self.help_resetDCS()

    def do_setRelayOn(self, arg):
        try:
            args = arg.split()
            self.socketCall( 'setRelayOn ' + args[0] )
        except:
            self.help_setRelayOn()

    def do_setRelayOff(self, arg):
        try:
            args = arg.split()
            self.socketCall( 'setRelayOff ' + args[0] )
        except:
            self.help_setRelayOff()

    def do_restoreInterlock( self, arg ):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            self.socketCall( 'restoreInterlock' )
        except:
            self.help_restoreInterlock()
        
    def do_restartPID( self, arg ):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            msg = self.socketCall( 'restartPID' )
            return msg
        except:
            self.help_restartPID()
        
    def do_stopPID( self, arg ):
        self.debugMsg( sys._getframe().f_code.co_name, arg )
        
        try:
            msg = self.socketCall( 'stopPID' )
            return msg
        except:
            self.help_stopPID()
            
    def do_setLV_On( self, arg ):
        while True:
            print( "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" )
            print( "WARNING: you are going to power-up the module!" )
            print( "If you are 100% sure about what you are going to do, please press 'y', otherwise any key pressing will candel the command." )
            print( "For example, the IV must be taken before LV is made ON." )
            print( "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" )
            ans = input( ">>> " ).lower()

            if ans == 'y':
                
                try:
                    msg = self.socketCall( 'setLV_On cli' )
                    return msg
                except:
                    self.help_setLV_On()

                break
            else:
                print( "LV is untouched." )
                break

        
    def do_setLV_Off( self, arg ):
        try:
            msg = self.socketCall( 'setLV_Off cli' )
            return msg
        except:
            self.help_setLV_Off()
        
    # def do_setHV_On( self, arg ):
    #     try:
    #         msg = self.socketCall( 'setHV_On cli' )
    #         return msg
    #     except:
    #         self.help_setHV_On()
        
    def do_setHV_Off( self, arg ):
        try:
            msg = self.socketCall( 'setHV_Off cli' )
            return msg
        except:
            self.help_setHV_Off()
        
    def do_resetMPOD( self, arg ):
        try:
            msg = self.socketCall( 'resetMPOD cli' )
            return msg
        except:
            self.help_resetMPOD()
        
    def do_setVerbosity(self, arg):
        try:
            args = arg.split()
            self.socketCall( 'setVerbosity ' + args[0] )
        except:
            self.help_setVerbosity()


    def emptyline(self):
        pass

    def help_kill_server(self):
        print ('This command kills the server, and closes the client console.')
        
    def help_relaunch_server(self):
        print ('This command kills the server, and relaunch a new server.\nOptions for pidServer.py can be appended (e.g. relaunch_server -t 20 -V -1.0)')
        
    def help_quit(self):
        print ('With Ctrl-D or \'quit\' or \'exit\', one can disconnect from the server, and exit the client console.')
        
    def help_EOF(self):
        print ('With Ctrl-D or \'quit\' or \'exit\', one can disconnect from the server, and exit the client console.')
        
    def help_exit(self):
        print ('With Ctrl-D or \'quit\' or \'exit\', one can disconnect from the server, and exit the client console.')
        
    def help_log(self):
        print ('log [n=1]: show the last n records of the monitor log.')
        
    def help_params(self):
        print ('params: show the current PID control parameters.')
        
    def help_setT(self):
        print ('setT <value>: set the target temperature to the specified value. Special enum are appointed')
        print ('  - exchange     : +30 degC')
        print ('  - warm         : +20 degC')
        print ('  - cold         : -15 degC')
        print ('  - cold_startup : -35 degC')
        
    def help_resetV(self):
        print ('resetV <value>: force reset the Peltier voltage to the specified value.')
        
    def help_setKp(self):
        print ('setKp <value>: change the coefficient of P-feedback to the specified value.')
        
    def help_setKi(self):
        print ('setKp <value>: change the coefficient of I-feedback to the specified value.')
        
    def help_setKd(self):
        print ('setKd <value>: change the coefficient of D-feedback to the specified value.')
        
    def help_setMaxP(self):
        print ('setMaxP <value>: change the ceiling of the P-feedback strength to the specified value.')
        
    def help_setMaxI(self):
        print ('setMaxI <value>: change the ceiling of the I-feedback strength to the specified value.')
        
    def help_setMaxD(self):
        print ('setMaxD <value>: change the ceiling of the D-feedback strength to the specified value.')
        
    def help_setIduration(self):
        print ('setIduration <value>: change the length of integral to the specified value.')
        
    def help_setTimeOffsetP(self):
        print ('setTimeOffsetP <value>: change the PID extrapolation time offset (positive = advance).')
        
    def help_setTimeOffsetI(self):
        print ('setTimeOffsetI <value>: change the PID extrapolation time offset (positive = advance).')
        
    def help_setTimeOffsetD(self):
        print ('setTimeOffsetD <value>: change the PID extrapolation time offset (positive = advance).')
        
    def help_restoreInterlock(self):
        print ('restoreInterlock: when interlock is active and the system state is out of the interlock dangerous range, by this command the interlock flag can be restored. Note that PID feedback restarting is only possible after restoration of the interlock flag.')
        
    def help_restartPID(self):
        print ('restartPID: when PID feedback is disabled, restart it manually.\nNote that PID does not activate automatically once stopped.\nRestarting will only succeed when the interlock is inactive.')
        
    def help_stopPID(self):
        print ('stopPID: when PID feedback is disabled, restart it manually. Note that PID does not activate automatically once stopped.')
        
    def help_unlock(self):
        print ('unlock: unlock the box cover and disable the interlock.')
        
    def help_resetDCS(self):
        print ('resetDCS: restore the Arduino hardware interlock status')
        
    def help_setRelayOn(self):
        print ('setRelayOn [1-8]: switch the relay to HIGH.')
        print ('  1: Heater')
        print ('  2: Peltier interlock')
        print ('  3: Chiller interlock')
        print ('  4: Peltier "+" line polarity (DO NOT USE USUALLY)')
        print ('  5: Peltier "-" line polarity (DO NOT USE USUALLY)')
        print ('  6: LV interlock')
        print ('  7: HV interlock')
        print ('  8: Box Lock Solenoid (DO NOT USE USUALLY)')
        
    def help_setRelayOff(self):
        print ('setRelayOff [1-8]: switch the relay to LOW')
        print ('  1: Heater')
        print ('  2: Peltier interlock')
        print ('  3: Chiller interlock')
        print ('  4: Peltier "+" line polarity (DO NOT USE USUALLY)')
        print ('  5: Peltier "-" line polarity (DO NOT USE USUALLY)')
        print ('  6: LV interlock')
        print ('  7: HV interlock')
        print ('  8: Box Lock Solenoid (DO NOT USE USUALLY)')
        
    def help_setLV_On(self):
        print( 'Set module\'s LV power channel to ON state' )
        
    def help_setLV_Off(self):
        print( 'Set module\'s LV power channel to OFF state' )
        
    # def help_setHV_On(self):
    #     print( 'Set module\'s HV power channel to ON state' )
        
    def help_setHV_Off(self):
        print( 'Set module\'s HV power channel to OFF state' )
        
    def help_resetMPOD(self):
        print( 'Clear MPOD event states like interlock signals' )
        
    def help_setVerbosity(self):
        print ('setVerbosity [info, debug, verbose]: switch the server log verbosity')
        
