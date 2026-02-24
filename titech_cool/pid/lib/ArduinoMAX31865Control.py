#/usr/bin/python3

import serial
import time

class ArduinoMAX31865Control:
    m_port = '/dev/ttyUSB0'
    m_baud = 9600

    def __init__(self, config):
        if 'port' in config:
            self.m_port = config['port']
        if 'baud' in config:
            self.m_baud = config['baud']
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        time.sleep(1)

    def getR(self):
        while True:
            self.m_ctrl.write( b'r' )
            line = self.m_ctrl.readline()
            try:
                return float( line )
            except:
                continue

    def getT(self):
        r = self.getR()
        T = -329. + 0.3593 * r + 2.9e-5 * r * r
        return T

if __name__ == '__main__':
    rtd = ArduinoMAX31865Control( { "port" : "/dev/ttyUSB3", "baud" : 9600 } )
    print( 'R = ', rtd.getR() )
    print( 'T = ', rtd.getT() )

