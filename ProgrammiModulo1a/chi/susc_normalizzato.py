import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pylab

L_list=[10,15,20,25,30,35,40,45,50]
data={L: pylab.loadtxt('chiL'+str(L)+'.txt',unpack = True) for L in L_list}

for L in L_list:
    data[L][0]=(data[L][0]-0.44)*L
    data[L][1]=data[L][1]*L**(-7/4)
    data[L][2]=data[L][2]*L**(-7/4)

colors = matplotlib.cm.jet(np.linspace(0,1,len(L_list)))
for index, L in enumerate(L_list):
    pylab.errorbar(data[L][0],data[L][1],data[L][2],linestyle = '', color = colors[index], marker = '.', label='L='+str(L))



pylab.rc('font')
pylab.xlabel('($\u03b2 - \u03b2 _c$)L')
pylab.ylabel('$\chi / L^{7/4}$')
pylab.title('Profilo della funzione $\phi$')
pylab.minorticks_on()
plt.grid()
plt.legend(loc='upper right')
pylab.show()
