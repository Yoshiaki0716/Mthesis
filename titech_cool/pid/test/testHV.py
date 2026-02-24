#!/usr/bin/env python3

import sys
import optparse
import traceback

from PidClient import *
from PidMonitor import *
from Logger import *

parser = optparse.OptionParser()
parser.add_option( '-c', '--command', action = 'store', default = 'status', choices=['start', 'stop', 'status'])
parser.add_option( '-v', '--verbose', action = 'store', default = 'info', dest = 'verbosity' )
parser.add_option( '--hostname', action = 'store', default = '127.0.0.1', dest = 'host' )
parser.add_option( '-p', '--port', action = 'store', type = 'int', default = 50007, dest = 'port' )
options, remainder = parser.parse_args( sys.argv )

logger = Logger( '{}/logs/client.log'.format( os.getenv('PIDPATH') ), verbosity = options.verbosity )
mon = PidMonitor(500)
cli = PidClient( logger, mon, options )

#--------------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        if options.command == 'start':
            msg = cli.socketCall( 'setLV_On cli' )
            msg = cli.socketCall( 'setHV_On cli' )
            cli.do_quit('')
        elif options.command == 'stop':
            msg = cli.socketCall( 'setHV_Off cli' )
            msg = cli.socketCall( 'setLV_Off cli' )
            cli.do_quit('')
        print('HV ON duration for the module: {}'.format('hoge'))
    except:
        print(traceback.format_exc())
        pass
