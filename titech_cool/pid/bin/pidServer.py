#!/usr/bin/env python3

import sys
import os
import traceback
import psutil
import datetime
import time
import daemon
import pprint

#--------------------------------------------------
# Find if any other 'pidServer.py' is running on the machine
# and the script only permits to run at most single process
#--------------------------------------------------

for proc in psutil.process_iter():
    try:
        processName = proc.name()
        processID = proc.pid
        if processName.find('pidServer') >= 0 and proc.pid != os.getpid():
            print( processName , ' ::: ', processID, ' is already running --> exit.' )
            sys.exit(1)
    except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
        pass


#--------------------------------------------------
# Multi-threading
#--------------------------------------------------
import threading

#--------------------------------------------------
# Option parsers
#--------------------------------------------------
import optparse

#--------------------------------------------------
# Json stuffs
#--------------------------------------------------
import json

#--------------------------------------------------
# Logger
#--------------------------------------------------
from Logger import *

#--------------------------------------------------
# Deploy low-level device controllers by config
#--------------------------------------------------
from DeviceInvoker import *

#--------------------------------------------------
# PID control and monitoring
#--------------------------------------------------
from PidControl import *
from PidMonitor import *


#--------------------------------------------------------------------------------

def configure(options, config, logger):
    
    #--------------------------------------------------
    # reading configs
    #--------------------------------------------------
    
    logger.info( 'pidServer init configuration options:')
    opt_dict = vars( options )
    for k in opt_dict:
        logger.info( '[config] {}: {}'.format( k, opt_dict[k] ) )
    #for k in config['pidControl']:
       # logger.info( '[config] {}: {}'.format( k, config['pidControl'][k] ) )
    
    logger.info('test1')
    pid = PidControl( config['pidControl'] )
    logger.info('test2')
    mon = PidMonitor( config['pidMonitor']['maxHistory'] )
    logger.info('test3')
    
    mon.logOutFile   = options.log
    pid.logOutFile   = options.log

    logger.info('test4')

    pid.T_set = options.target
    pid.T_set_regulated = options.target

    logger.info('test5')

    #############
    #logger.info('%s', pid)
    ############
    #--------------------------------------------------
    # instantiating device controllers
    #--------------------------------------------------
    logger.info('test5.5')
    
    devices = invokeDevices( config['devices'] )
    
    logger.info('test6')

    peltCtrl = devices['peltierControl']

    logger.info('test7')

    dcsCtrl  = devices['dcsControl']
    
    peltCtrl.setOn()
    peltCtrl.setOff()
    peltCtrl.resetAlarms()
    
    logger.info('test8')

    peltCtrl.setOn()
    time.sleep(3)
    try:
        if abs( options.Vinit ) < 10.0 :
            if options.Vinit < 0.0:
                peltCtrl.setV( 0.0 )
                time.sleep(2)
                dcsCtrl.pelPolInverted()
            peltCtrl.setV( abs(options.Vinit) )
            pid.V     = options.Vinit
            pid.Vprev = options.Vinit
    except:
        print(traceback.format_exc())
        logger.fatal( 'could not communicate with the peltier control power supply device. terminating' )
        print( 'could not communicate with the peltier control power supply device. terminating' )
        exit(1)
        
    time.sleep(3)
    pid.injectDependencies( mon, devices, logger, config )
    
    return pid, devices
    
    logger.info('test9')
#--------------------------------------------------------------------------------
def main_unit():

    now = datetime.datetime.now()
    now_unix = int(time.mktime(now.timetuple()))
    today = now.date().isoformat()
    wd = os.getenv('PIDPATH')

    with open( f'{wd}/config/default.json') as f:
        defaults = json.load(f)
    
    #--------------------------------------------------
    parser = optparse.OptionParser()
    log_dir = os.environ["EQC_OPERATOR_DIR"]+"/../log/"+os.environ["HOSTNAME"]
    os.makedirs( log_dir, exist_ok = True )
    
    log_path = log_dir+"/"+os.environ["HOSTNAME"]+"_"+today+"_log.txt"

    parser.add_option( '-c', action = 'store', default = os.path.join( defaults['config'] ), dest = 'configPath' )
    parser.add_option( '-o', action = 'store', default = log_path, dest = 'log' )
    parser.add_option( '-V', action = 'store', default = -9999, dest = 'Vinit', type = 'float' )
    parser.add_option( '-t', '--target', action = 'store', default = 20.0, dest = 'target', type = 'float'  )
    parser.add_option( '-v', '--verbose', action = 'store', default = 'info', dest = 'verbosity' )
    options, remainder = parser.parse_args( sys.argv )
    options.configPath = wd + '/' + options.configPath
    
    with open( options.configPath ) as f:
        config = json.load(f)
    
    logger = Logger( options.log, verbosity = options.verbosity, use_dest = False )
    
    logger.info( '==================================================' )
    logger.info( 'service launched' )
    logger.info( '==================================================' )
    #####
    #logger.info('test')
    ####
    
    #--------------------------------------------------
    # instantiations and configurations
    #--------------------------------------------------

    pid, devices = configure(options, config, logger)
    #logger.debug( 'Finished configuration.')
    logger.info( 'Finished configuration.')
    #####
    if config["pidControl"]["LVtype"] == "CUSTOM":
        logger.info( 'test10' )
    #####
    if config["pidControl"]["LVtype"] == "RIGOL":
        os.system( "killall -9 watchdog_rigol.sh" )
        os.system( "${PIDPATH}/bin/watchdog_rigol.sh >/var/tmp/watchdog_rigol.log &" )
        logger.info( 'RIGOL watchdog restarted.')
    
    #--------------------------------------------------
    # start process threads
    #--------------------------------------------------

    tCommand = threading.Thread( target = pid.listenCommands, args = (50007,) )
    tCommand.start()
    logger.info( 'Finished thread.')

    logger.info( 'started to wait for client connection' )

            
    #--------------------------------------------------
    # Exception handling in the following loop
    #
    # - 'kill' flag is used for user's intention of ending the service
    #    ==> then the threads are joined and break the loop
    #
    # - 'fatal' flag is used for system's unexpected exception raising
    #    ==> then the threads end, but the master service process
    #        attempts to re-launch all threads and continue running.
    #--------------------------------------------------
    
    while not pid.fatal:
        
        try:

            if not pid.kill:

                tPID = threading.Thread( target = pid.pidLoop )
                tPID.start()
                logger.info( 'started PID loop thread' )
                logger.info( 'service is ready' )

                pid.sendEmail( 'Started PID loop', 'Started PID loop as the server initialization.\nService is in operation', category='INFO' )
                
                #--------------------------------------------------
                # end process
                #--------------------------------------------------
                
                tCommand.join()
                tPID.join()

                if pid.kill:
                    break
                    
        except:
            
            pid.fatal = True
            tPID.join()
            logger.error( 'fatal error was detected.' )
            logger.error( 'restarting device control threads' )
            pid.fatal = False

            
    #--------------------------------------------------
    # Termination processes
    #--------------------------------------------------
    
    # devices['dcsControl'].pelPolNormal()
    # devices['peltierControl'].resetAlarms()
    # devices['peltierControl'].setV(0)
    # devices['peltierControl'].setOff()
    
    logger.info( 'Service is properly terminating.' )
    logger.info( '==================================================' )
    logger.info( '' )
    logger.info( '' )
    logger.info( '' )
    exit(0)


if __name__ == '__main__':
    
    wd = os.getenv('PIDPATH')
    now = datetime.datetime.now()
    now_unix = int(time.mktime(now.timetuple()))
    
    with open( f'{wd}/config/default.json' ) as f:
        defaults = json.load(f)
    
    parser = optparse.OptionParser()
    parser.add_option( '-c', action = 'store', default = os.path.join( wd, defaults['config'] ), dest = 'configPath', help='specify the control config file of the service.' )
    parser.add_option( '-o', action = 'store', default = '{}/logs/{}.log'.format( wd, now_unix ), dest = 'log', help='use this option if you like to use a custom log output' )
    parser.add_option( '-V', action = 'store', default = -9999, dest = 'Vinit', type = 'float', help = 'set the initial Peltier voltage in [V]' )
    parser.add_option( '-t', '--target', action = 'store', default = 20.0, dest = 'target', type = 'float', help = 'set the initial target set temperature [degC]'   )
    parser.add_option( '-v', '--verbose', action = 'store', default = 'info', dest = 'verbosity', help = 'set the initial verbosity of the server. Options are { info, debug, verbose }. You can change it later by the client.' )
    options, remainder = parser.parse_args( sys.argv )
    options.configPath = wd + '/' + options.configPath
    
    dc = daemon.DaemonContext(
        stderr=open('/tmp/pidServer_console.txt', 'w+'),
        stdout=open('/tmp/pidServer_console.txt', 'w+')
    )
    
    with dc:
        main_unit()
