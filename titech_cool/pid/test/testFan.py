#!/usr/bin/python

import RPi.GPIO as GPIO
import threading
import time

pin = 3
fraction = 0.1
duration = 0.02
kill = False

GPIO.setwarnings(False)
GPIO.setmode(GPIO.BOARD)
GPIO.setup(pin, GPIO.OUT)

def fanControl():
    global kill
    while True:
        GPIO.output(pin, True)
        time.sleep(duration * fraction)
        GPIO.output(pin, False)
        time.sleep(duration * (1.0 - fraction) )
        
        if kill == True:
            print "killed"
            GPIO.output(pin, False)
            return
            

if __name__ == '__main__':
    t = threading.Thread( target = fanControl )
    t.start()

    while True:
        s = raw_input('> ')
        if s == 'kill' :
            print "kill is called"
            kill = True
            time.sleep(0.1)
            break
        else:
            if( float( s ) < 0 ): continue
            if( float( s ) > 1 ): continue
            fraction = float( s )
    
    t.join()
