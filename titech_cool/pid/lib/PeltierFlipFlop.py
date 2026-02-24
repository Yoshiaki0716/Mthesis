#!/usr/bin/python

import RPi.GPIO as GPIO
import time

class PeltierFlipFlop:
    def __init__(self, config):
        self.pin1 = config['pins'][0]
        self.pin2 = config['pins'][1]
        self.status = False
        
        GPIO.setwarnings(False)
        GPIO.setmode(GPIO.BOARD)
        GPIO.setup(self.pin1, GPIO.OUT)
        GPIO.setup(self.pin2, GPIO.OUT)
        GPIO.output(self.pin1, False)
        GPIO.output(self.pin2, False)
    '''
    def __init__(self, pin1 = 7, pin2 = 10):
        self.pin1 = pin1
        self.pin2 = pin2
        self.status = False
        
        GPIO.setwarnings(False)
        GPIO.setmode(GPIO.BOARD)
        GPIO.setup(self.pin1, GPIO.OUT)
        GPIO.setup(self.pin2, GPIO.OUT)
        GPIO.output(self.pin1, False)
        GPIO.output(self.pin2, False)
    '''
        
    def flip(self):
        self.status = not self.status
        GPIO.output(self.pin1, self.status)
        time.sleep( 0.1 )
        GPIO.output(self.pin2, self.status)

    def read(self):
        return [ GPIO.input(self.pin1), GPIO.input(self.pin2) ]
        
    def setNormal(self):
        self.status = False
        GPIO.output(self.pin1, self.status)
        time.sleep( 0.1 )
        GPIO.output(self.pin2, self.status)
        
    def setInverted(self):
        self.status = True
        GPIO.output(self.pin1, self.status)
        time.sleep( 0.1 )
        GPIO.output(self.pin2, self.status)
        

