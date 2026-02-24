#/usr/bin/python3

import RPi.GPIO as GPIO
import time
import sys

class MAX31855:

    spi_ss=3
    spi_clk=19
    spi_miso=21
    spi_mosi=23
    wait=0.01
    
    def __init__(self,ss,clk=19,miso=21,mosi=23,w=0.01):

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


        GPIO.output(self.spi_clk, False)
        time.sleep(self.wait)
        
        #initialize
        GPIO.output(self.spi_ss,   False)
        time.sleep(self.wait)

        sign=1
        out=self.read()
        if(out==1):
            sign=-1

        temp=0
        for k in range(0,13):
            temp<<=1
            out=self.read()
            if(out==1):
                temp|=0x1

        #14 reserved
        out = self.read()

        #15 Fault
        out = self.read()
        #if(out==1):
        #    print "Error at 15"

        #16 sign for ref
        out = self.read()

        
        signR=1
        if(out==1):
            signR=-1

        tempR=0
        for k in range(0,11):
            tempR<<=1
            out = self.read()
            if(out==1):
                tempR |= 0x1

        #28 reserved
        out = self.read()

        #29
        out = self.read()
        #if(out==1):
        #    print "Error : Short TC to V"
        #30
        out = self.read()
        #if(out==1):
        #    print "Error : Short TC to G"
            
        #31
        out = self.read()
        #if(out==1):
        #    print "Error : Open TC"
            
        #end process
        GPIO.output(self.spi_clk, False)
        time.sleep(self.wait)
        GPIO.output(self.spi_ss,   True)
        time.sleep(self.wait)
        GPIO.output(self.spi_clk, True)
        time.sleep(self.wait)
        if(sign==1):
            return sign*temp*0.25
        else:
            return sign*(2**11-temp*0.25)




class ThermocoupleCalibConfig:
    def __init__( self, params = [ 0.0006892, -0.05584, 1.99 ] ):
        self.params = params

    def getT( self, adc ) :
        return adc - (self.params[0]*pow(adc,2) + self.params[1] * adc + self.params[2])
        

class ThermocoupleReader:
    def __init__( self, config ):
        self.device  = MAX31855( ss = config['pinout']['ss'] )
        self.channel = config['channel']
        self.calib   = ThermocoupleCalibConfig( config['calib'] )
        
    def getT(self):
        return self.calib.getT( self.device.readValue( self.channel ) )


if __name__ == '__main__':
    thermocouple1 = ThermocoupleReader( { 'pinout' : { 'ss' : 13 }, 'channel' : 0, 'calib': [ 0.0006892, -0.05584, 1.99 ] } )
    thermocouple2 = ThermocoupleReader( { 'pinout' : { 'ss' : 15 }, 'channel' : 0, 'calib': [ 0.0006892, -0.05584, 1.99 ] } )
    print( thermocouple1.getT() )
    print( thermocouple2.getT() )
