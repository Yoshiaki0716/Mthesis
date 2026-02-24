from datetime import datetime
import os
import time

class Logger:
    levels = { 'fatal': 0, 'error': 1, 'warning': 2, 'info': 3, 'debug': 4, 'verbose': 5 }
    def __init__( self, dest, verbosity = 'info', use_dest = True ):
        self.dest = dest
        self.verbosity = self.levels[verbosity]
        self.use_dest = use_dest
        
    def write( self, verbosity, message ):

        now = datetime.now()
        
        if self.use_dest:
            log_path = self.dest
        else:
            today = now.date().isoformat()
            log_dir = os.environ["EQC_OPERATOR_DIR"]+"/../log/"+os.environ["HOSTNAME"]
            log_path = log_dir+"/"+os.environ["HOSTNAME"]+"_"+today+"_log.txt"
            
        with open( log_path, "a" ) as f:
            f.write( ' '.join( [ str( int(time.time() ) ), now.strftime("%Y-%m-%d %H:%M:%S"), "[ pidServer ]", '{:10s}'.format(verbosity), message, '\n'] ) )

    def setVerbosity( self, verbosity ):
        self.verbosity = self.levels[verbosity]

    def fatal( self, message ):
        if self.verbosity >= self.levels['fatal']: self.write( 'FATAL', message )
        
    def error( self, message ):
        if self.verbosity >= self.levels['error']: self.write( 'ERROR', message )
        
    def warning( self, message ):
        if self.verbosity >= self.levels['warning']: self.write( 'WARNING', message )
        
    def info( self, message ):
        if self.verbosity >= self.levels['info']: self.write( 'INFO', message )
        
    def debug( self, message ):
        if self.verbosity >= self.levels['debug']: self.write( 'DEBUG', message )
        
    def verbose( self, message ):
        if self.verbosity >= self.levels['verbose']: self.write( 'VERBOSE', message )
