#/usr/bin/python3

import serial
import time
import math

class ArduinoNTCReader:
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
            self.m_ctrl.write( b'ModuleTemp\n' )
            line = self.m_ctrl.readline().decode().rstrip()
            try:
                return float( line )
            except:
                return -100.0

    def getMUX(self):
        while True:
            self.m_ctrl.write( b's\n' )
            line = self.m_ctrl.readline().decode().rstrip()
            try:
                if line == "No Connection" or line == "":
                    return [float(-1000.0),float(-1000.0),float(-1000.0),float(-1000.0)]
                else:
                    values = line.split()
                    if len(values)==4:
                        values=[float(i) for i in values]
                        return values
                    else:
                        return [float(-1000.0),float(-1000.0),float(-1000.0),float(-1000.0)]
            except:
                return [float(-1000.0),float(-1000.0),float(-1000.0),float(-1000.0)]



if __name__ == '__main__':
    ntc = ArduinoNTCReader( { "port" : "/dev/ttyACM0", "baud" : 9600 } )
    print( 'temp = ', ntc.getT() )
    MUX = ntc.getMUX()
    print(MUX)
    for i in range(0,4):
        print( "MUX_channel"+ str(i)+ ": "+ str(MUX[i]))
        print( type(MUX[i]))
