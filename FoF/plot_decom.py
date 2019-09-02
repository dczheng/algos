#!/usr/bin/env python

from matplotlib import use
use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np
import sys

ds = []
N = int(sys.argv[1])
for i in range(N):
    fn = "p_%03i.dat"%i
    ds.append( np.loadtxt( fn ) )

xmin = ymin = 1e30
xmax = ymax = -1e30

fig, ax = plt.subplots( 1, 1, figsize=(8,8) )
for i in range(N):
    x = ds[i][:,0]
    y = ds[i][:,1]
    xmin = np.min( [xmin, x.min()] )
    ymin = np.min( [ymin, y.min()] )
    xmax = np.max( [xmax, x.max()] )
    ymax = np.max( [ymax, y.max()] )
    ax.plot( x, y )
print( xmin, xmax, ymin, ymax )

ax.set_xlim( [xmin, xmax] )
ax.set_ylim( [ymin, ymax] )

plt.savefig( 'decomp.png' )
