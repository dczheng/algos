#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

ds = []
N = 8
for i in range(N):
    fn = "p_%03i.dat"%i
    ds.append( np.loadtxt( fn ) )

if N%2:
    x = N//2+1
else:
    x = N//2
fig, axs = plt.subplots( x, 2, figsize=(4*2,4*x) )

xmin = ymin = 1e30
xmax = ymax = -1e30

for i in range(N):
    ax = axs[i//2, i%2]
    x = ds[i][:,0]
    y = ds[i][:,1]
    xmin = np.min( [xmin, x.min()] )
    ymin = np.min( [ymin, y.min()] )
    xmax = np.max( [xmax, x.max()] )
    ymax = np.max( [ymax, y.max()] )
    ax.plot( x, y )
print( xmin, xmax, ymin, ymax )

for i in range(N):
    ax = axs[i//2, i%2]
    ax.set_xlim( [xmin, xmax] )
    ax.set_ylim( [ymin, ymax] )

plt.savefig( 'decomp.png' )
