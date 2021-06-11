
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pylab
L_list=[10,15,20,25,30,35,40,45,50]
data={L: pylab.loadtxt('eneL'+str(L)+'.txt',unpack = True) for L in L_list}

colors = matplotlib.cm.jet(np.linspace(0,1,len(L_list)))
for index, L in enumerate(L_list):
    pylab.errorbar(data[L][0],data[L][1],data[L][2],linestyle = '', color = colors[index], marker = '.', label='L='+str(L))


#plt.ylim(0, 1)
#plt.xlim(0.35, 0.51)
pylab.rc('font')
pylab.xlabel('$\u03b2$', fontsize=20)
pylab.ylabel('$\epsilon$', fontsize=20)
pylab.minorticks_on()
plt.grid()
plt.legend(loc='upper right')
pylab.show()
