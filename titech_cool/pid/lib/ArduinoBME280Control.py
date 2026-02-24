#/usr/bin/python3

import serial
import time

class ArduinoBME280Control:
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
            self.m_ctrl.write(b'T')
            line = self.m_ctrl.readline()
            try:
                return float( line )
            except:
                continue

    def getRH(self):
        while True:
            self.m_ctrl.write(b'H')
            line = self.m_ctrl.readline()
            try:
                return float( line )
            except:
                continue

    def getP(self):
        while True:
            self.m_ctrl.write(b'P')
            line = self.m_ctrl.readline()
            try:
                return float( line )
            except:
                continue



if __name__ == '__main__':
    bme280 = ArduinoBME280Control( { "port" : "/dev/ttyACM0", "baud" : 9600 } )
    print( 'temp = ', bme280.getT() )
    print( 'rh = ', bme280.getRH() )
    print( 'p = ', bme280.getP() )
