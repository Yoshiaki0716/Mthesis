#!/usr/bin/env python3

import sys
import optparse
import traceback
from cmd import Cmd

import socket
import threading
import time
import re
import json
from PidControl import *
from PidClient import *
from PidMonitor import *
from Logger import *

parser = optparse.OptionParser()
parser.add_option( '-v', '--verbose', action = 'store', default = 'info', dest = 'verbosity' )
parser.add_option( '--hostname', action = 'store', default = '127.0.0.1', dest = 'host' )
parser.add_option( '-p', '--port', action = 'store', type = 'int', default = 50007, dest = 'port' )
options, remainder = parser.parse_args( sys.argv )

logger = Logger( '{}/logs/client.log'.format( os.getenv('PIDPATH') ), verbosity = options.verbosity )
mon = PidMonitor(500)
cli = PidClient( logger, mon, options )

def getMon():
    global cli
    global mon
    global logger
    
    while True:

        if cli.kill: break
        
        try:
            
            while True:
                if cli.isBusy:
                    print( 'waiting for other thread unlocks the socket...' )
                    time.sleep(0.5)
                else:
                    break
            
            cli.isBusy = True
            
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect( (options.host, options.port) )
            
            command = 'sync all' if len( mon.data ) == 0 else ('sync ' + str( mon.data['elapsed'][-1] ) )
            logger.debug( 'sending command ' + command )
            
            sock.sendall( command.encode() )
            
            logger.debug( 'sync command sent...' )
            
            buf = b''
            while True:
                data = sock.recv(1024)
                buf += data

                logger.verbose( 'received data size: ' + str( len(data) ) )
                
                if repr(data).find( '\\n' ) < 0:
                    continue
                else:
                    break
            
            sock.close()
            cli.isBusy = False
            
            logger.debug( 'completed data reception. size = : ' + str( len(buf) ) )
            
            j = json.loads( buf.decode("utf-8") )

            if len(j) == 0:
                pass
            else:
                mon.t0 = j[0]
                for k, v in j[1].items():
                    if k in mon.data:
                        mon.data[k] += v
                    else:
                        mon.data.update( { k:v } )
                mon.status = j[2]

            for i in range(5):
                time.sleep(1)
                if cli.kill: break
        
            
        except Exception as e:

            if cli.kill: return
            
            # logger.error( 'getMon() Exception: ' + str(e) )
            # print( str(e) )
            # print(traceback.format_exc()) 
            time.sleep(5)
            
        cli.isBusy = False



#--------------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        
        #tMon = threading.Thread( target = getMon )
        #tMon.start()
        
        while True:
            
            cli.cmdloop()
            
            if cli.kill == True:
                break
            
        #tMon.join()
        
    except:
        print(traceback.format_exc())
        pass
    
