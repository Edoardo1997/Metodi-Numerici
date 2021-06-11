import matplotlib
import matplotlib.pyplot as plt
import pylab
import numpy
from scipy.optimize import curve_fit

#x1,y1,dy1=pylab.loadtxt('binder7k_0.35beta.txt',unpack=True)
x2, y2, dy2 = pylab.loadtxt('Binder_vs_L_beta0.50.txt',unpack = True)
#pylab.errorbar(x1,y1,dy1,linestyle = '', color = 'navy', marker = '.', label='L=20')
pylab.errorbar(x2,y2,dy2,linestyle = '', color = 'red', marker = '.')



# bellu

# AT THE FIRST ATTEMPT COMMENT FROM HERE TO THE END

# define the function (linear, in this example)
def ff(x, aa, bb,cc):
    return aa+bb*x**2+cc*x**4

# define the initial values (STRICTLY NEEDED!!!)
init=(3,-102,2260)

# prepare a dummy xx array (with 2000 linearly spaced points)




# plot the best fit curve
pylab.grid()
#pylab.plot(xx,ff(xx,*pars), color='red')
pylab.rc('font',size=18)
pylab.xlabel('$L$', fontsize=20)
pylab.ylabel('$B$', fontsize=20)
pylab.title('Cumulante di Binder per $\u03b2=0.50$')
pylab.minorticks_on()
#plt.legend(loc='upper right')
pylab.show()
