#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

d = np.loadtxt( 'hidx2.dat' )
x = d[:,1]
y = d[:,2]
print( x.min(), x.max() )
print( y.min(), y.max() )
fig, ax = plt.subplots( 1, 1, figsize=(6,6) )
ax.plot( x, y )
fig.savefig( "hidx.png" )
plt.close()

d = np.loadtxt( 'hidx3.dat' )
x = d[:,1]
y = d[:,2]
z = d[:,3]
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot( x, y, z )
plt.show()

