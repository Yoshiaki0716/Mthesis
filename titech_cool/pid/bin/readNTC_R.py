#!/usr/bin/env python3
import time
import sys
import optparse
import numpy as np
import scipy.optimize as op
import pandas as pd
import math
from scipy.optimize import fsolve
SteinHart_A = 0.8676453371787721e-3
SteinHart_B = 2.541035850140508e-4
SteinHart_C = 1.868520310774293e-7

from InfluxAccess import *

ctrl = InfluxAccess( { "host" : '192.168.100.104', "port" : 8086, "database" : "REPIC",
                       "measurement" : "CoolingBox01",
                       "field" : "Temperature (ch.Module) [C]" } )

Temp = ctrl.getLast()
Temp_abs = Temp + 273.15
#print( ctrl.getLast() )



def f(x):
    return SteinHart_A + SteinHart_B*math.log(x) + SteinHart_C*math.log(x)*math.log(x)*math.log(x)

def g(x):
    return 1/Temp_abs

def h(x):
    return f(x) - g(x)

results = fsolve(h,1)
Regi = results[0]
print(Regi)
#print("results : ",results)
#print("R : ",Regi)
#print("T [C] : ",ctrl.getLast() )
#print("f(x) : ",f(results))
#print("g(x) : ",g(results))
