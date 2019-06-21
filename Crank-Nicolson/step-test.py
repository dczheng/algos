#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use( 'ggplot' )


x = np.logspace( -3, 0, 100 )

y = np.exp(x) - 1
y1 = 1 - np.exp(-x)
y2 = 1 - np.exp(-2*x)

plt.plot( x, y, label=r'$e^{\Delta}-1$' )
plt.plot( x, y1,'--', label=r'$1-e^{-\Delta}$' )
plt.plot( x, y2, '-.', label=r'$1-e^{-2\Delta}$' )
plt.legend()
plt.xlabel( r'$\Delta$' )
plt.xscale( 'log' )

plt.show()

