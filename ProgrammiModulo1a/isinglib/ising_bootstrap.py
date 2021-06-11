#!/usr/bin/env python3

'''Bootstrap algorithm, Metropolis stories and observables calculation'''

from isinglib import ising_type
import typing as tp
import numpy as np


rng = np.random.default_rng()


def bootstrap(osservabile: ising_type.tpobs, vec: tp.List[float], boot_cycle: int = 100) -> float:
    '''Bootstrap function with increasing binning'''
    bootk = []
    len_vec = len(vec)
    #Starts with bin > 1, to reach the asymptotic behaviour faster
    bin_vec = 1 + len_vec//1000

    #Max number of iterations: 6
    while(bin_vec <= 1 + len_vec/10):
        Obs=[]
        #Resampling
        for _ in range(boot_cycle):
            temp = []
            for _ in range(int(len_vec/bin_vec)):
                i = rng.integers(len_vec)
                temp.extend(vec[i:min(i+bin_vec,len_vec)])
                
            Obs.append(osservabile(temp))
        bootk.append(np.std(Obs))

        #Increasing binning exponentially
        bin_vec*=2
    return max(bootk)



def media_abs(vec: tp.List[float], valass: bool = False) -> float:
    '''Calculates mean value of a vector, with absolute value of elements if required'''
    if valass:
        return np.mean([abs(i) for i in vec])
    else:
        return np.mean(vec)


def varianza_abs(vec: tp.List[float], valass: bool = False) -> float:
    '''Calculates standard deviation of a vector, with absolute value of elements if valass=True'''
    if valass:
        return np.var([abs(i) for i in vec])
    else:
        return np.var(vec)


def binder(vec: tp.List[float]) -> float:
    '''Calculates binder cumulant'''
    m2=0
    m4=0
    for i in vec:
        m2+=i**2
        m4+=i**4
    return len(vec)*m4/m2**2



def punto(vec: tp.Dict[str, tp.List[float]], L: int ,  name: str ) -> tp.Dict[str, float] :
    '''Calculates mag (magnetizzation), amag (absolute magnetization value), chi (magnetic susceptibility), e (energy), c (specific heat) and binder (binder cumulant). 
    Energy and magnetization are the mean values of the corresponding stories
    Specific heat is proportional to the variance of energy story 
    Susceptibility is proportional to the variance of the absolute value of magnetizations in the magnetization story
    Returns the resulting point (value and error of the observable "name") in the dictionary res'''
    
    res = {}
    quant = name in ['amag','chi','mag','binder']
    #Quant = True if the observable requires magnetization stories, False if requires energy stories. 

    if(name in ['amag','ene','mag']):
        #mag and ene do not require absolute value, amag requires absolute value of magnetizations
        if(name == 'mag'):
            quant = False 
        res['valore'] = media_abs(vec, quant) 
        res['errore'] = bootstrap(lambda x: media_abs(x, quant), vec) #error calculation

    #susc requires absolute value, c does not require absolute value
    elif(name in ['c','chi']):
        res['valore'] = L*L*varianza_abs(vec, quant) 
        res['errore'] = bootstrap( lambda x: L*L*varianza_abs(x, quant), vec)

    elif(name == 'binder'):
        res['valore'] = binder(vec)
        res['errore'] = bootstrap(binder, vec)

    return res

