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



#Enter in second mode, plot only a previous simulation
L_list=[10,15,20,25,30,35,40,45,50]
for L in L_list:
    out_file = f'binder_chi_amag_mag_c_ene_L{L}_h0.000.txt'
    datas = salva.read_data(out_file)
    for oss in datas['d_oss']:
        if oss not in os.listdir(os.curdir):
            os.mkdir(oss)
        f=open(oss+os.sep+oss+f"L{datas['L']}"+'.txt','w')
        for x, y, dy in zip(datas['x_axis'], datas['d_oss'][oss]['valore'], datas['d_oss'][oss]['errore']):
            if 0.34<=x<0.481:
                print(f'{x:8f} {y:8f} {dy:8f}', file=f)
        f.close()

