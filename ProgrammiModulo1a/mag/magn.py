
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pylab
L_list=[10,20,30,40,50,60]
data={L: pylab.loadtxt('magL'+str(L)+'.txt',unpack = True) for L in L_list}

colors = matplotlib.cm.jet(np.linspace(0,1,len(L_list)))
for index, L in enumerate(L_list):
    pylab.errorbar(data[L][0],data[L][1],data[L][2],linestyle = '', color = colors[index], marker = '.', label='L='+str(L))




pylab.rc('font')
pylab.xlabel('$\u03b2$', fontsize=20)
pylab.ylabel('$M$', fontsize=20)
pylab.minorticks_on()
plt.grid()
plt.legend(loc='upper right')
pylab.show()
