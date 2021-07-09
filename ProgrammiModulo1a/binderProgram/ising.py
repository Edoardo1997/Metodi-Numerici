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
    nstep=30000
    i_dec=50
    iflag=1

    L_list= list(range(210,310,10))
   # for x in [0.45, 0.35, 0.50]:
    for x in [0.50]:
        bind = {'valore': [], 'errore': []}
        if f'Binder_beta{x:.2f}' not in os.listdir(os.curdir):
            os.mkdir(f'Binder_beta{x:.2f}')
            print('New Directory for Binder_vs_L created')
        for L in L_list:
            seed = np.random.randint(1000000000)
            #save input for C program
            iflag = int(x<=0.44)
            f=open('input.txt','w')
            print(f"{seed}\n{L}\n{x:.4f}\n{extfield:.3f}\n{nstep}\n{i_dec}\n{iflag}", file=f)
            f.close()
            fmt = f'binderL{L}'+os.extsep+'txt'
    
            #out_file = "{}_L{}_h{:.3f}".format('_'.join(supp_obs),L,extfield) + os.extsep + 'txt'
            
            res = sp.run('./story ', shell=True) #Run C program to calculate stories
            
 
            #Dictionary of stories
            v = { 'ene' : [] , 'magn' : [] }
        
            f=open('misure.txt','r')
            line=f.readline()
            while line.strip():
                s = line.strip()
                try:
                    v['magn'].append(float(s))
                except (ValueError, TypeError, IndexError):
                    print(f'Incorrect measures at the output of C program at temperature {x:4f} with L={L}')
                    sys.exit(1)
                line=f.readline()
            f.close()

            salva.salva_storia(L, x, extfield , vec = v, file_name = os.curdir + os.sep + f'Binder_beta{x:.2f}' + os.sep +fmt)

           
            P = bts.punto(v['magn'], L, name = 'binder')
            print(x, L, P['valore'], P['errore'], sep=' ')
            bind['valore'].append(P['valore'])
            bind['errore'].append(P['errore'])
        
        f = open(f'Binder_vs_L_beta{x:.2f}.txt', 'w')
        for L, b, db in zip(L_list, bind['valore'], bind['errore']):
            print(f'{L} {b:.16f} {db:.16f}', file=f)
        f.close()
        

       

    




