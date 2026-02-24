import math
import serial
import time

class TakasagoKXControl:
    def __init__(self, config):
        self.m_port = config['port']
        self.m_addr = config['address']
        if 'baud' in config:
            self.m_baud = config['baud']
        else:
            self.m_baud = 9600

        time.sleep(1)
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
            
    def setOn(self):
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'A{}\n'.format( self.m_addr ).encode() )
        self.m_ctrl.write( 'OT1\n'.encode() )
        time.sleep(0.5)

    def setOff(self):
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'A{}\n'.format( self.m_addr ).encode() )
        self.m_ctrl.write( 'OT0\n'.encode() )
        time.sleep(0.5)

    def setV(self, outV):
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)

        #prev = self.getV().get('meas')
        prevobj = self.getV()
        prev =  prevobj.get('meas')
        prev_set = prevobj.get('set')
        I = self.getI()
        Ilim = I.get('limit')

        print(f"takasagotest:lim={Ilim},Vmeas={prev},Vset={prev_set},set={outV}")

        currentLimit = outV * 0.8
        
        if currentLimit < 0.7 : currentLimit = 0.7
        if currentLimit > 15.40 : currentLimit = 15.40

        if outV > prev:
            print("case1")
            self.m_ctrl.write( 'A1,LC{:.2f},LV{:.2f}\n'.format( currentLimit, outV*1.2 ).encode() )
            self.m_ctrl.write( 'A1,OC{:.2f},OV{:.2f}\n'.format( currentLimit*0.9, outV ).encode() )
        else:
            print("case2")
            self.m_ctrl.write( 'OC{:.2f},OV{:.2f}\n'.format( currentLimit*0.9, outV ).encode() )
            self.m_ctrl.write( 'LC{:.2f},LV{:.2f}\n'.format( currentLimit, outV*1.2 ).encode() )
        time.sleep(0.05)
    def getStatus(self):
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'TK4\n'.encode() )
        line = self.m_ctrl.readline().decode()
        if line.find('STAT') < 0:
            raise ValueError( 'Read Error in TakasagoKXControl.getStatus() : ' + line )
        if line.find( 'ALM128' ) >= 0:
            raise ValueError( line )
        
        return { 'P_ON'  :int(line[10]),
                 'OHP'   :int(line[9]),
                 'LIMIT' :int(line[8]),
                 'OCP'   :int(line[7]),
                 'OVP'   :int(line[6]),
                 'CC'    :int(line[5]),
                 'CV'    :int(line[4])   }

    def getSetting(self):
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'TK0\n'.encode() )
        line = self.m_ctrl.readline().decode()
        
        if line.find( 'ALM128' ) >= 0:
            raise ValueError( line )
        if len(line) == 0:
            raise ValueError( 'null string passed' )

        return { 'SetVoltage'   :float(line.split(',')[0]),
                 'SetCurrent'   :float(line.split(',')[1]),
                 'VoltageLimit' :float(line.split(',')[2]),
                 'CurrentLimit' :float(line.split(',')[3]),
                 'Output'       :int(line.split(',')[4]),
                 'SINK'         :int(line.split(',')[5])    }

    def resetAlarms(self):
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'AR1\n'.encode() )
        time.sleep(0.5)
                       
    def getV(self):
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'TK5\n'.encode() )
        time.sleep(0.05)
        line = self.m_ctrl.readline().decode()
        if line == '' :
            raise ValueError( 'Read Error in TakasagoKXControl.getV()' )
        if line.find( 'ALM128' ) >= 0:
            raise ValueError( line )
        measValue = float( line.split(',')[0].replace('V','') )

        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'TK0\n'.encode() )
        line = self.m_ctrl.readline().decode()
        
        if line.find( 'ALM128' ) >= 0:
            raise ValueError( line )
        if len(line) == 0:
            raise ValueError( 'null string passed' )
        setValue = float(line.split(',')[0])
        limit = float(line.split(',')[2])
        return { 'set': setValue, 'meas': measValue, 'limit': limit }

    def getI(self):
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'TK5\n'.encode() )
        line = self.m_ctrl.readline().decode()
        measValue = float( line.split(',')[1].replace('A', '') )

        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'TK0\n'.encode() )
        line = self.m_ctrl.readline().decode()
        if line.find( 'ALM128' ) >= 0:
            raise ValueError( line )
        if len(line) == 0:
            raise ValueError( 'null string passed' )
        setValue = float(line.split(',')[1])
        limit = float(line.split(',')[3])
        return { 'set': setValue, 'meas': measValue, 'limit': limit }


class TakasagoZXControl:
    m_port = '/dev/ttyUSB0'
    m_addr = 1
    m_baud = 9600

    def __init__(self, config):
        self.m_port = config['port']
        self.m_addr = config['address']
        if 'baud' in config:
            self.m_baud = config['baud']
        else:
            self.m_baud = 9600
        
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)
        self.m_ctrl.write( 'System:COMMUNICATE:SERIAL:PACE OFF\n'.encode() )
        self.m_ctrl.readline()

    def setOn(self):
        self.m_ctrl.write( 'ADDR {}\n'.format( self.m_addr ).encode() )
        self.m_ctrl.write('OUTPUT ON\n'.encode() )

    def setOff(self):

        self.m_ctrl.write( 'ADDR {}\n'.format( self.m_addr ).encode() )
        self.m_ctrl.write('OUTPUT OFF\n'.encode() )

    def setV(self, outV):
        n=2    
        voltage = math.floor(outV * 10 ** n)/(10 ** n)
        self.m_ctrl.write( 'VOLT {:.2f}\n'.format( outV ).encode() )


    def getID(self):
        self.m_ctrl.write( '*idn?\n'.encode() )
        msg = self.m_ctrl.readline().decode()[:-1]
        return msg
        
    def getStatus(self):
        self.m_ctrl.write( 'STATUS:MEASURE:CONDITION?\n'.encode() )
        try:
            out = int( self.m_ctrl.readline().decode()[:-1] )
            
            status = { 'EXT_TRIP_LT' : out & (1<<17),
                       'EXT_TRIP'    : out & (1<<16),
                       'OVP'         : out & (1<<15),
                       'OCP'         : out & (1<<14),
                       'CP'          : out & (1<<13),
                       'EXT_ON'      : out & (1<<12),
                       'SYS_ALM'     : out & (1<<11),
                       'DD_ON'       : out & (1<<10),
                       'MST/BST'     : out & (1<<9),
                       'P-ON_B'      : out & (1<<8),
                       'P-ON_M'      : out & (1<<7),
                       'AD-OV_ALM'   : out & (1<<6),
                       'OHP_ALM'     : out & (1<<5),
                       'OHP_ALM'     : out & (1<<4),
                       'OVP_ALM'     : out & (1<<3),
                       'PL'          : out & (1<<2),
                       'CC'          : out & (1<<1),
                       'CV'          : out & (1<<0)        }
            
            
            return status
        except:
            return {}

    def resetAlarms(self):
        self.m_ctrl.write( 'ALM:CLEar\n'.encode() )
        time.sleep(0.5)
                       
    def getV(self):
        self.m_ctrl.write( 'MEASURE:VOLTAGE?\n'.encode() )
        line = self.m_ctrl.readline().decode()
        if line == '' :
            raise ValueError( 'Read Error in TakasagoZXControl.getV()' )
        if line.find( 'ERROR' ) >= 0:
            raise ValueError( line )
        return float( line )

    def getI(self):
        self.m_ctrl.write( 'MEASURE:CURRENT?\n'.encode() )
        line = self.m_ctrl.readline().decode()
        return float( line )


class TexioPFRControl:
    def __init__(self, config):
        self.m_port = config['port']
        self.m_addr = config['address']
        if 'baud' in config:
            self.m_baud = config['baud']
        else:
            self.m_baud = 9600

        time.sleep(1)
        self.m_ctrl = serial.Serial(self.m_port, self.m_baud, timeout=1)

    def idn(self):
        try:
            self.m_ctrl.write( '*idn?\n'.encode() )
            line = self.m_ctrl.readline().decode().rstrip()
            print( line )
        except:
            raise Exception('Communication error with TexioPFR Power Supply')
            
    def setOn(self):
        try:
            self.m_ctrl.write( 'OUTP ON\n'.encode() )
        except:
            raise Exception('Communication error with TexioPFR Power Supply')

    def setOff(self):
        try:
            self.m_ctrl.write( 'OUTP OFF\n'.encode() )
        except:
            raise Exception('Communication error with TexioPFR Power Supply')

    def setV(self, value, limit = None):
        if limit == None:
            limit = min( 12.0, value*1.2 )
            
        value = min( 10.0, value )
            
        try:
            self.command( f'volt:prot {limit}', None )
            self.command( f'volt {value}', None )
        except:
            raise Exception('Communication error with TexioPFR Power Supply')
        
        time.sleep(0.1)
        
        return self.getV()

    def setI(self, value, limit = None):
        limit = min( value*1.2, 10.0 )
        if limit == None:
            limit = min( 10.0, value*1.2 )
            
        value = min( 9.0, value )
            
        try:
            self.command( f'curr:prot {limit}', None )
            self.command( f'curr {value}', None )
        except:
            raise Exception('Communication error with TexioPFR Power Supply')
        
        time.sleep(0.1)
        
        return self.getI()

    def getStatus(self):
        try:
            return self.command( f'output?', bool )
        except:
            raise Exception('Communication error with TexioPFR Power Supply')
        
    def isTripped(self):
        # print( self.command( 'output:protection:tripped?', bool ) )
        Vinfo = self.getV()
        Iinfo = self.getI()
        
        is_acceptable = Iinfo.get('meas', 0) < 4.0 if Vinfo.get('meas') < 1.1 else abs( Iinfo.get('meas') - 1.00 * (Vinfo.get('meas',0) - 1.47) ) < 4.0
        
        return not is_acceptable

    def getSetting(self):
        try:
            self.m_ctrl.write( 'APPL?\n'.encode() )
            line = self.m_ctrl.readline().decode().rstrip()
            return line
        except:
            raise Exception('Communication error with TexioPFR Power Supply')

    def resetAlarms(self):
        try:
            self.m_ctrl.write( ':OUTP:PROT:CLE\n'.encode() )        
        except:
            raise Exception('Communication error with TexioPFR Power Supply')
                       
    def getV(self):
        try:
            setValue  = self.command( 'volt?', float )
            measValue = self.command( 'meas:volt?', float )
            limit     = self.command( 'volt:prot?', float )
        except:
            raise Exception('Communication error with TexioPFR Power Supply')
            
        return { 'set': setValue, 'meas': measValue, 'limit': limit }

    def getI(self):
        try:
            setValue  = self.command( 'curr?', float )
            measValue = self.command( 'meas:curr?', float )
            limit     = self.command( 'curr:prot?', float )
        except:
            raise Exception('Communication error with TexioPFR Power Supply')
            
        return { 'set': setValue, 'meas': measValue, 'limit': limit }
        
    def command(self, command, outType = str ):
        self.m_ctrl.write( f'{command}\n'.encode() )
        if command.find("?")>0:
            line = self.m_ctrl.readline().decode().rstrip()
            return outType( line )
        else:
            return None
        
