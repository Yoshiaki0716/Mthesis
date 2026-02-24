import RPi.GPIO as GPIO
import time
import sys
import numpy as np

class MCP3208:

    def __init__(self,ss = 7, clk = 19, miso =21, mosi = 23, w = 0.01):

        self.spi_ss=ss
        self.spi_clk=clk
        self.spi_miso=miso
        self.spi_mosi=mosi
        self.wait=w

        GPIO.setwarnings(False)
        GPIO.setmode(GPIO.BOARD)
        GPIO.setup(self.spi_mosi, GPIO.OUT)
        GPIO.setup(self.spi_miso, GPIO.IN)
        GPIO.setup(self.spi_clk, GPIO.OUT)
        GPIO.setup(self.spi_ss,  GPIO.OUT)

        GPIO.output(self.spi_ss,   False)
        time.sleep(self.wait)

        GPIO.output(self.spi_ss,   True)
        time.sleep(self.wait)

        GPIO.output(self.spi_clk, True)
        time.sleep(self.wait)

    def write(self,val):
        GPIO.output(self.spi_clk,  False)
        time.sleep(self.wait)
        if(val==0):
            GPIO.output(self.spi_mosi, False)
        else:
            GPIO.output(self.spi_mosi, True)
        time.sleep(self.wait)
        GPIO.output(self.spi_clk,  True)
        time.sleep(self.wait)

    def read(self):
        GPIO.output(self.spi_clk, False)
        time.sleep(self.wait)
        GPIO.output(self.spi_clk, True)
        time.sleep(self.wait)
        out=0
        if (GPIO.input(self.spi_miso)):
            out=1
        return out

    def readValue(self,ch):
        #initialize
        GPIO.output(self.spi_ss,   False)
        time.sleep(self.wait)

        #start bit,single ended
        self.write(1)
        self.write(1)

        #select channel
        assert (ch >= 0 and ch < 8), 'Error : no such channel: {0}'.format( ch )
        
        for bit in range(3):
            self.write( ( ch & (1 << bit) ) >> bit )

        #wait one clock
        self.write(1)
			

        out=self.read()
        if(out==1):
            print( "Error: should be null" )
		
        value = 0
        for k in range(0,12):
            value <<= 1
            out=self.read()
            if (out==1):
                value |= 0x1
        #end process
        GPIO.output(self.spi_clk, False)
        time.sleep(self.wait)
        GPIO.output(self.spi_ss,   True)
        time.sleep(self.wait)
        GPIO.output(self.spi_clk, True)
        time.sleep(self.wait)

        return value



class NtcCalibConfig:
    def __init__( self, R25 = 100., Rext = 470., B = 4281., T0 = (273.15 + 27.26) ):
        self.R25 = R25
        self.Rext = Rext
        self.B = B
        self.T0 = T0

    def getT( self, adc ) :
        r = adc / 4096.
        R = self.Rext * r / ( 1. - r )
        T = 1. / ( np.log( R / self.R25 ) / self.B + 1. /self.T0 )
        return  T - 273.15
        

class NtcReader:
    def __init__( self, config ):
        pinout = config['pinout']
        calibs = config['calib']
        self.device  = MCP3208( ss = pinout['ss'], clk = pinout['clk'], miso = pinout['miso'], mosi = pinout['mosi'] )
        self.channel = config['channel']
        self.calib   = NtcCalibConfig( R25 = calibs['R25'], Rext = calibs['Rext'], B = calibs['B'], T0 = calibs['T0'] )
        
    '''
    def __init__( self, device, channel, calib ):
        self.device = device
        self.channel = channel
        self.calib = calib
    '''
        
    def getT(self):
        adc = self.device.readValue( self.channel )
        return self.calib.getT( adc )
