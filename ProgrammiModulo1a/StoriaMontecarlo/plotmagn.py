
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pylab

x1, y1= pylab.loadtxt('magn_step_20L_Beta0.35.txt',unpack = True)
x2, y2= pylab.loadtxt('magn_step_30L_Beta0.35.txt',unpack = True)
x3, y3= pylab.loadtxt('magn_step_40L_Beta0.35.txt',unpack = True)
x4, y4= pylab.loadtxt('magn_step_50L_Beta0.45.txt',unpack = True)

xu1, yu1= pylab.loadtxt('magn_step_20L_Beta0.35u.txt',unpack = True)
xu2, yu2= pylab.loadtxt('magn_step_30L_Beta0.35u.txt',unpack = True)
xu3, yu3= pylab.loadtxt('magn_step_40L_Beta0.35u.txt',unpack = True)
xu4, yu4= pylab.loadtxt('magn_step_50L_Beta0.45u.txt',unpack = True)

xd1, yd1= pylab.loadtxt('magn_step_20L_Beta0.35d.txt',unpack = True)
xd2, yd2= pylab.loadtxt('magn_step_30L_Beta0.35d.txt',unpack = True)
xd3, yd3= pylab.loadtxt('magn_step_40L_Beta0.35d.txt',unpack = True)
xd4, yd4= pylab.loadtxt('magn_step_50L_Beta0.45d.txt',unpack = True)


pylab.errorbar(x4,y4,linestyle = '', color = 'blue', marker = '.')
pylab.errorbar(xd4,yd4,linestyle = '', color = 'red', marker = '.')
pylab.errorbar(xu4,yu4,linestyle = '', color = 'green', marker = '.')
plt.legend(loc='upper right')
plt.xlabel('Tempo monte carlo')
plt.ylim(-1, 1)


# pylab.errorbar(x1,y1,linestyle = '', color = 'navy', marker = '|', label='L=20')
# pylab.errorbar(x2,y2,linestyle = '', color = 'orchid', marker = '|', label='L=30')
# pylab.errorbar(x3,y3,linestyle = '', color = 'aquamarine', marker = '|', label='L=40')
# pylab.errorbar(x4,y4,linestyle = '', color = 'darkgreen', marker = '|', label='L=50')
#
#
#
# plt.title('Magnetizzazione media per \u03b2 = 0.35')
# plt.xlabel('Tempo monte carlo')
plt.ylabel('Magnetizzazione')

#plt.legend(loc='upper right')
pylab.show()
