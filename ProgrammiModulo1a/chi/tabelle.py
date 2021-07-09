
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pylab
from scipy.optimize import curve_fit

    

f = open("suscettivitamax.txt", "r")
line=f.readline()
line=f.readline()
while line:
    vec=[float(i) for i in line.split()]
    print('{vec[0]:.0f} & {vec[1]:.4f} $\pm$ {vec[2]:.4f} & {vec[3]:.2f} $\pm$ {vec[4]:.2f} & {vec[5]:.2f}\\\\'.format(vec=vec))
    line=f.readline()



f.close()






