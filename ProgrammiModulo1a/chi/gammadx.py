
import matplotlib
import matplotlib.pyplot as plt
import pylab
import numpy
from scipy.optimize import curve_fit

x,y,Dy=pylab.loadtxt('chiL50.txt',unpack=True)

x1,y1,Dy1=pylab.loadtxt('chiL50.txt',unpack=True)
x=[]
y=[]
Dy=[]
bc=0.4324
for i in range(len(x1)):
    if 0.4395<= float(x1[i]) <= 0.4515:
        x.append(x1[i])
        y.append(y1[i])
        Dy.append(Dy1[i])
x=numpy.array(x)
y=numpy.array(y)
Dy=numpy.array(Dy)
x=-bc/x+1
x=numpy.log(x)
Dy=Dy/y
y=numpy.log(y)
print(x)

pylab.errorbar(x,y,Dy,linestyle = '', color = 'black', marker = '.')



z=0


# bellurie
pylab.rc('font',size=18)
pylab.xlabel('$log(|t|)$')
pylab.ylabel('$log(\chi)$')
pylab.title('Suscettività magnetica su reticolo 50x50 in scala logaritmica')
pylab.minorticks_on()

xx=numpy.linspace(min(x),max(x),2000)
# AT THE FIRST ATTEMPT COMMENT FROM HERE TO THE END

# define the function (linear, in this example)

def ff(x, aa, bb):
    return aa-bb*x

# define the initial values (STRICTLY NEEDED!!!)
init=(-2.5,-1.75)

# prepare a dummy xx array (with 2000 linearly spaced points)


# plot the fitting curve computed with initial values
# AT THE SECOND ATTEMPT THE FOLLOWING LINE MUST BE COMMENTED
#pylab.plot(xx,ff(xx,*init), color='blue')

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
print('gamma/nu=', pars[1],'con incertezza',numpy.sqrt(covm[1,1]))

# plot the best fit curve
pylab.plot(xx,ff(xx,*pars), color='red')

# show the plot
pylab.show()
