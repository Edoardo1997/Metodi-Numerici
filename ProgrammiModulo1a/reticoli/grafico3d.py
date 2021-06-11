import matplotlib
import matplotlib.pyplot as plt
import pylab
import numpy as np
from scipy.optimize import curve_fit

L=50
z035=pylab.loadtxt('reticolo035.txt',unpack=True)
z040=pylab.loadtxt('reticolo040.txt',unpack=True)
z043=pylab.loadtxt('reticolo043.txt',unpack=True)
z045=pylab.loadtxt('reticolo045.txt',unpack=True)
i=[int(x/L) for x in range(0,L*L)]
j=[x%L for x in range(0,L*L)]





plt.figure(figsize=(6, 4))

#subplot(2, 2, 1)
a1=plt.subplot(2, 2, 1,projection = '3d')
a1.scatter3D(i, j, z035)
# Give labels
a1.set_xlabel('i')
a1.set_ylabel('j')
a1.set_zlabel('spin')
a1.set_title('$\u03b2=0.35$')

a2=plt.subplot(2, 2, 2,projection = '3d')
a2.scatter3D(i, j, z040)
# Give labels
a2.set_xlabel('i')
a2.set_ylabel('j')
a2.set_zlabel('spin')
a2.set_title('$\u03b2=0.40$')


a3=plt.subplot(2, 2, 3,projection = '3d')
a3.scatter3D(i, j, z043)
# Give labels
a3.set_xlabel('i')
a3.set_ylabel('j')
a3.set_zlabel('spin')
a3.set_title('$\u03b2=0.43$')


a4=plt.subplot(2, 2, 4,projection = '3d')
a4.scatter3D(i, j, z045)
# Give labels
a4.set_xlabel('i')
a4.set_ylabel('j')
a4.set_zlabel('spin')
a4.set_title('$\u03b2=0.45$')

plt.show()
