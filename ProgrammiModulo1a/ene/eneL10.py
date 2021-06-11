
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pylab
x0, y0, dy0 = pylab.loadtxt('eneL10allpoints.txt',unpack = True)

x1 = np.linspace(0, 0.20, 1000)
y1 = -2*x1
plt.title('L=10    h=0.000', fontsize=20)
pylab.errorbar(x1,y1, color = 'red', marker = '', label='f($\u03b2$)=-2$\u03b2$')
pylab.errorbar(x0,y0,dy0,linestyle = '', color = 'black', marker = '.', label='$\epsilon$')






pylab.rc('font')
pylab.xlabel('$\u03b2$',fontsize=20)
pylab.ylabel('$\epsilon$',fontsize=20)
pylab.minorticks_on()
plt.grid()
plt.legend(loc='upper right')
pylab.show()
