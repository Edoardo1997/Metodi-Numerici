import numpy as np
import matplotlib.pyplot as plt
import random as rand
from matplotlib import cm

# sample data

"""
columns=3
data=1000

A=np.zeros((columns,data))

ranges = 5

limits=[]
for i in range(ranges-1):
    limits.append(rand.randint(1,data-1))


s_limit = sorted(limits)
s_limit.append(data)


K = np.arange(ranges)+1

start = 0
for index,limit in enumerate(s_limit):
    print(limit)
    for i in range(limit-start):
        A[0,start+i]=K[index]
    start=limit

A[1,:] = np.arange(data)
A[2,:] = A[]
"""

ones, N, g, h, m, N_max, mx, mz= np.loadtxt('/Users/edoardo/Metodi_Numerici/results.txt',unpack=True)

def prepare_plotter(color_index, x, y):
    colors = np.unique(color_index)
    #print(colors)
    curves=[]
    for color in colors:
        i_color=int(color)
        #print(i_color)
        #print(colors==i_color)
        x_curve = x[color_index==i_color]
        y_curve = y[color_index==i_color]
        curves.append([x_curve, y_curve])
    return curves, colors

curves, color_curve = prepare_plotter(N, g, mz)
np_curves = np.asarray(curves)

x_curves = np_curves[:,0,:]
y_curves = np_curves[:,1,:]


fig, ax = plt.subplots()  # Plotting data
ax.set_title('Magnetization along z axis')
ax.set_xlabel('$g$')
ax.set_ylabel('$M_z$')


low_lim = 15
big_lim = 25


low_lim-=1
n = big_lim-low_lim
colors = cm.jet(np.linspace(0,1,n))


for index,mz in enumerate(y_curves):
    if(index<low_lim or index>=big_lim):
        continue
    ax.plot(x_curves[index, :], mz, color=colors[index-low_lim],label='N = {}'.format(int(color_curve[index]))) 
    ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax.legend(fontsize=7)
plt.show()

    