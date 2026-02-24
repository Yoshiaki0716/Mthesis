import RPi.GPIO as GPIO
import time


class CpuFanControl:
    pwmRate = 0.1
    pwmPeriod = 0.02

    def __init__(self, config):
        self.pin = config['pin']
        
        GPIO.setwarnings(False)
        GPIO.setmode(GPIO.BOARD)
        GPIO.setup(self.pin, GPIO.OUT)
        GPIO.output(self.pin, False)
        
    def bind( self, pid ):
        self.pid = pid
        
    def fanControl(self):
        
        while True:
            GPIO.output(self.pin, True)
            time.sleep(self.pwmPeriod * self.pwmRate)
            GPIO.output(self.pin, False)
            time.sleep(self.pwmPeriod * (1.0 - self.pwmRate) )
            
            if self.pid.kill == True:
                GPIO.output(self.pin, False)
                # print "CpuFanControl.fanControl(): killed"
                return
            
