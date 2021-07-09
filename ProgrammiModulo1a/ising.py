#!/usr/bin/env python3

'''2D ising model, main program
Calculates mag (magnetizzation), amag (absolute magnetization value), chi (magnetic susceptibility), e (energy), c (specific heat) and binder (binder cumulant)
from a simulation, using Metropolis algotithm
'''

from isinglib import ising_bootstrap as bts
from isinglib import ising_files as salva
from isinglib import ising_small as sml
from isinglib import ising_plot as grf
import subprocess as sp
import numpy as np
import sys
import os


#Check if the user has the correct python version
version = sys.version_info.major >= 3 and sys.version_info.minor >= 6
if not version:
    print('Python 3.6 or later version is needed.')
    sys.exit(3)

#List of supported observables:
supp_obs = ['binder','chi','amag','mag','c','ene']


print('{:#^61}\n\n{}'.format('  classical ising 2D simulation project  '.title(), 'Choose program mode:'))

    #user choice between mode
while True:
    temp = input('1: simulation\t\t\t2: plot the previous results\n').lower().strip()
    if temp in ['1','','2'] :
        break
    elif temp == 'q':
        print('Quitting')
        sys.exit(0)
    else:
        print('Not understood, try again. Press q to quit')

mod = (temp in ['1',''])

#Enter in second mode, plot only a previous simulation
if not mod:
    out_file = input('You are currently in:\n{}\nInsert file path/filename or filename if it is in the cwd \n'.format(os.getcwd())).strip()
    try:
        if input('Choose 1 or 2\n1: all plots in one figures\t\t\t2: one figure for each plot\n').strip()=='1':
            grf.mode_3(out_file)
        else:
            grf.mode_2(out_file)
    except (FileNotFoundError, IsADirectoryError) as e:
        print('Error found while plotting previous results: ')
        print(str(e).split(']')[1].strip())
        sys.exit(1)

    #Raise if there is unexpected or missing information
    except (IndexError, ValueError):
        print(f"File {out_file} does not match standard results files format")
        sys.exit(2)

#Simulation mode
else:
    extfield=0
    nstep=10000
    i_dec=100
    #create x axis
    x_axis=[]
    beta = 0.34
    while beta < 0.481:
        x_axis.append(beta)
        if 0.40 <= beta < 0.45 :
            beta+=0.0005
        else:
            beta+=0.001
            
            
    if 'MC_stories' not in os.listdir(os.curdir):
        os.mkdir('MC_stories')
        print('Directory MC_stories created')
    for L in range(15,52,5):
        path= os.curdir
        out_file = "{}_L{}_h{:.3f}".format('_'.join(supp_obs),L,extfield) + os.extsep + 'txt'
        
       
        if f"L={L}" not in os.listdir(os.curdir + os.sep + 'MC_stories'):
            os.mkdir(os.curdir+os.sep+'MC_stories'+os.sep+f"L={L}" )
            print(f"Directory L={L} created" )
        print('MonteCarlo stories will be saved into the directory ''{}L={}'''.format(os.curdir+os.sep+'MC_stories'+os.sep ,L))
        #Creates list of beta with given grain and observables dictionary
      
        d_oss = { oss : {'valore': [], 'errore': []} for oss in supp_obs }


        iflag=1 #configuration all random
        #Calculates values corresponding to given beta range
        for x in x_axis:
            seed = np.random.randint(1000000000)
            #save input for C program
            f=open('input.txt','w')
            print(f"{seed}\n{L}\n{x:.4f}\n{extfield:.3f}\n{nstep}\n{i_dec}\n{iflag}", file=f)
            f.close()
            iflag=3 #no new lattice is created because the previous one (at the previous value of temperature) is almost termalized at the current value of temperature (x)
            
           
            
            #filename
            fmt = 'h{h:.3f}_beta{beta:.4f}'.format(h=extfield,beta=x).replace('.',',')+os.extsep+'txt'
            


            res = sp.run('./story ', shell=True) #Run C program to calculate stories
            
       

            #Dictionary of stories
            v = { 'ene' : [] , 'magn' : [] }
        
            f=open('misure.txt','r')
            line=f.readline()
            while line.strip():
                s = line.split()
                try:
                    v['magn'].append(float(s[0]))
                    v['ene'].append(float(s[1]))
                except (ValueError, TypeError, IndexError):
                    print(f'Incorrect measures at the output of C program at temperature {x:4f} with L={L}')
                    sys.exit(1)
                line=f.readline()
            f.close()

            salva.salva_storia(L, x, extfield , vec = v, file_name = os.curdir + os.sep + 'MC_stories' + os.sep + f"L={L}"+os.sep+fmt)

            #Calculates value and error at the end of the simulation for every requested observable
            for oss in supp_obs:
                P = bts.punto(v[sml.find_matter(oss)], L, name = oss)
                d_oss[oss]['valore'].append(P['valore'])
                d_oss[oss]['errore'].append(P['errore'])
                if oss=='chi':
                    print(L,x,P['valore'],P['errore'])
        salva.func_save(L, extfield, 'beta', x_axis, d_oss, out_file, path)

        #If the user asked for only one value of beta (or T), prints the results on stdout, otherwise a plot is generated

        #for oss in supp_obs:
         #   block = oss == supp_obs[-1] #blocking last plot
          #  grf.plot_graph(x_axis, d_oss[oss]['valore'], d_oss[oss]['errore'], L ,extfield, 'beta', oss, block_fig = block)


    




