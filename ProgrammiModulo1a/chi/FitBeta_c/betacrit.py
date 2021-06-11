
import matplotlib
import matplotlib.pyplot as plt
import pylab
import numpy
from scipy.optimize import curve_fit

x,y,Dy,a,b=pylab.loadtxt('suscettivitamax.txt',unpack=True)


pylab.errorbar(x,y,Dy,linestyle = '', color = 'black', marker = '.')




# bellurie
pylab.rc('font',size=18)
pylab.xlabel('$L$  ')
pylab.ylabel('cavolo  ')
pylab.minorticks_on()

# AT THE FIRST ATTEMPT COMMENT FROM HERE TO THE END

# define the function (linear, in this example)
def ff(x, aa, bb, cc):
    return aa+bb*(x**(-1/cc))

# define the initial values (STRICTLY NEEDED!!!)
init=(0.44,-0.4,1)

# prepare a dummy xx array (with 2000 linearly spaced points)
xx=numpy.linspace(min(x),max(x),2000)

# plot the fitting curve computed with initial values
# AT THE SECOND ATTEMPT THE FOLLOWING LINE MUST BE COMMENTED
'''pylab.plot(xx,ff(xx,*init), color='blue')'''

# set the error
sigma=Dy
w=1/sigma**2

# call the minimization routine
pars,covm=curve_fit(ff,x,y,init,sigma, absolute_sigma=False)

# calculate the chisquare for the best-fit function
chi2 = ((w*(y-ff(x,*pars))**2)).sum()

# determine the ndof
ndof=len(x)-len(init)


# print results on the console
print('pars:',pars)
print('covm:',covm)
print ('chi2, ndof:',chi2, ndof)


print('beta_critico=', pars[0],'con incertezza',numpy.sqrt(covm[0,0]))
print('xbar=', pars[1],'con incertezza',numpy.sqrt(covm[1,1]))
print('nu=', pars[2],'con incertezza',numpy.sqrt(covm[2,2]))
# plot the best fit curve
pylab.plot(xx,ff(xx,*pars), color='red')


# show the plot
pylab.show()
