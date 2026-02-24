#!/usr/bin/python

import RPi.GPIO as GPIO
import time

pin = 7
pin2 = 10

GPIO.setwarnings(False)
GPIO.setmode(GPIO.BOARD)
GPIO.setup(pin, GPIO.OUT)
GPIO.setup(pin2, GPIO.OUT)

flag = True

def switch(period, tempo = 1.2, pin = 7):
    global flag
    if flag == False:
        GPIO.output(pin, True)
        flag = True
    else:
        GPIO.output(pin, False)
        flag = False
    
    time.sleep(period * tempo)

tempo = 1

while True:
    switch( 0.25 )
    switch( 0.125 )
    switch( 0.125 )
    switch( 0.125/3. )
    switch( 0.125/3. )
    switch( 0.125/3. )
    switch( 0.125 )
    switch( 0.125 )
    switch( 0.125 )
    switch( 0.5 )

    switch( 0.125, pin = 10 )
    switch( 0.125, pin = 10 )
    switch( 0.25, pin = 10 )
    switch( 0.25, pin = 10 )
    switch( 0.125/3., pin = 10 )
    switch( 0.125/3., pin = 10 )
    switch( 0.125/3., pin = 10 )
    switch( 0.125, pin = 10 )
    switch( 0.5, pin = 10 )

    switch( 0.25 )
    switch( 0.25 )
    switch( 0.25 )
    switch( 0.125 )
    switch( 0.125 )
    switch( 0.25 )
    switch( 0.25 )

    switch( 0.125, pin = 10 )
    switch( 0.125, pin = 10 )
    switch( 0.25, pin = 10 )
    switch( 0.25, pin = 10 )
    switch( 0.125/3., pin = 10 )
    switch( 0.125/3., pin = 10 )
    switch( 0.125/3., pin = 10 )
    switch( 0.125, pin = 10 )
    switch( 0.5, pin = 10 )

