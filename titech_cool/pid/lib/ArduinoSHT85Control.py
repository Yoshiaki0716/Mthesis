#/usr/bin/python3

import serial
import time
import math

class ArduinoSHT85Control:
    m_port = '/dev/ttyACM0'
    m_baud = 9600

    def __init__(self, config):
        if 'port' in config:
            self.m_port = config['port']
        if 'baud' in config:
            self.m_baud = config['baud']
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        time.sleep(1)

    def getT(self):
        while True:
            self.m_ctrl.write( b't' )
            line = self.m_ctrl.readline()
            try:
                return float( line )
            except:
                continue

    def getRH(self):
        while True:
            self.m_ctrl.write( b'h' )
            line = self.m_ctrl.readline()
            try:
                return float( line )
            except:
                continue

    def getDP_primitive(self):
        T = self.getT()
        RH = self.getRH()

        # Simple Magnus formula
        b = 18.678
        c = 257.14
        d = 234.5
        
        g = math.log(RH/100.) + b*T/(c+T)
        dp = c * g / ( b - g )
        
        return dp


    def getDP_standard(self):
        T = self.getT()
        RH = self.getRH()

        a = 5.64
        b = 0.861 * T - 46.4
        dp = a * math.sqrt(RH) + b
        
        return dp

    
    def getDP_conservative1sigma(self):
        T = self.getT()
        RH = self.getRH()

        a  = 0.00895 * math.exp( -0.317 * T ) + 7.60
        b  = 1.01 * T - 52.1 - 2.28e-9 * math.exp( -T )
        dp = a * math.sqrt(RH) + b
        return dp

    
    def getDP_conservative2sigma(self):
        T = self.getT()
        RH = self.getRH()

        a =  0.0135 * math.exp( -0.311*T ) + 8.32
        b =  0.832 * T -47.6 + 1.32e-10 * math.exp( -T )
        dp = a * math.sqrt(RH) + b

        return dp

    
    def getDP_conservative3sigma(self):
        T = self.getT()
        RH = self.getRH()
        
        a =  0.0183 * math.exp( -0.307*T ) + 9.03
        b =  0.653*T - 43.0 + 2.55e-9 * math.exp(-T)
        dp = a * math.sqrt(RH) + b
        
        return dp


if __name__ == '__main__':
    sht85 = ArduinoSHT85Control( { "port" : "/dev/ttyACM2", "baud" : 9600 } )
    print( 'temp = ', sht85.getT() )
    print( 'rh = ', sht85.getRH() )
    print( 'dew (primitive) = ', sht85.getDP_primitive() )
    print( 'dew (standard) = ', sht85.getDP_standard() )
    print( 'dew (conservative 1sigma) = ', sht85.getDP_conservative1sigma() )
    print( 'dew (conservative 2sigma) = ', sht85.getDP_conservative2sigma() )
    print( 'dew (conservative 3sigma) = ', sht85.getDP_conservative3sigma() )
