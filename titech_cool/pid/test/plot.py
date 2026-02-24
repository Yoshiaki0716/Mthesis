#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import math

x = []
y1 = []
y2 = []

plt.ioff()

fig, ( (ax1, ax2), (ax3, ax4) ) = plt.subplots( 2, 2, figsize=(16,12) )
fig.suptitle('Test Plots')

lines1, = ax1.plot( x,  y1 )
lines2, = ax2.plot( x,  y2 )

while True:
    
    if len(x) > 0:
        x.append( x[-1] + 1.e-1 )
    else:
        x.append( 0 )
        
    y1.append( math.sin(x[-1]) )
    y2.append( math.cos(x[-1]) )
    
    lines1.set_data( x, y1 )
    lines2.set_data( x, y2 )
    minx = max( 0, max(x) - 2*math.pi )
    maxx = max( 2*math.pi, max(x) * 1.1 )
    miny1 = min(y1) - max( 0.1, (max(y1)-min(y1)) * 0.2 )
    maxy1 = max(y1) + max( 0.1, (max(y1)-min(y1)) * 0.2 )
    miny2 = min(y2) - max( 0.1, (max(y2)-min(y2)) * 0.2 )
    maxy2 = max(y2) + max( 0.1, (max(y2)-min(y2)) * 0.2 )
    ax1.set_xlim(minx, maxx)
    ax1.set_ylim(miny1, maxy1)
    ax2.set_xlim(minx, maxx)
    ax2.set_ylim(miny2, maxy2)
    plt.pause(0.5)



