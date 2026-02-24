#/usr/bin/python3

import serial
import time
from multiprocessing import Process, Value

class ArduinoFlipper:
    m_port = '/dev/ttyACM0'
    m_baud = 9600

    def __init__(self, config):
        if 'port' in config:
            self.m_port = config['port']
        if 'baud' in config:
            self.m_baud = config['baud']
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.readline()
        self.status = False
        self.lockStatus = Value('i', -1 )
        self.setNormal()

    def flip(self):
        self.status = not self.status
        if self.status:
            self.m_ctrl.write( b'ON\n' )
        else:
            self.m_ctrl.write( b'OFF\n' )
        self.m_ctrl.readline()

    def setNormal(self):
        self.status = False
        self.m_ctrl.write( b'OFF\n' )
        self.m_ctrl.readline()
        
    def setInverted(self):
        self.status = True
        self.m_ctrl.write( b'ON\n' )
        self.m_ctrl.readline()

    def unlock(self):
        self.m_ctrl.write( b'UNLOCK\n' )
        time.sleep(0.1)

    def monitorLock(self):
        self.m_ctrl.write( b'LOCKSTAT\n' )
        str = ''
        while True:
            str = self.m_ctrl.readline().decode()
            str = str.replace('\r', '').replace('\n', '')
            if( str != '' ): break

        self.lockStatus.value = int(str)

        
if __name__ == '__main__':
    flipper = ArduinoFlipper( { "port" : "/dev/ttyUSB0", "baud" : 9600 } )
    flipper.flip()
    flipper.flip()
    flipper.monitorLock()
    print( flipper.lockStatus.value )
    flipper.unlock()
    flipper.monitorLock()
    print( flipper.lockStatus.value )
