#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

d = np.loadtxt( 'map.dat' )
gs = open( 'groups.dat' ).readlines()

fig, axs = plt.subplots( 2, 1, figsize=(5, 2*5) )
axs[0].imshow( d )

dd = np.zeros( d.shape ) 
xs = []
ys = []
idxs = []
for g in gs:
    t = [ int(tt) for tt in g.split() ]
    idxs.append(t[0])
    t = t[1:]
    xx = 0
    yy = 0
    gl = len(t) // 2
    for i in range( gl ):
        dd[ t[i*2+1], t[i*2] ] = 1
        xx += t[i*2]
        yy += t[i*2+1]
    xx /= gl
    yy /= gl
    xs.append( xx )
    ys.append( yy )

axs[1].imshow( dd )
#for i in range( len(xs) ):
#    axs[1].text( xs[i], ys[i], '%i'%idxs[i] )

plt.savefig( 'fof.png' )
