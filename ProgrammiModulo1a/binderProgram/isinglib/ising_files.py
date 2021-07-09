#!/usr/bin/env python3

'''Module to handle files'''

from isinglib import ising_errors as errors
from isinglib import ising_small as sml
from isinglib import ising_type
import typing as tp
import os



def func_save(L: int, h: float, x_name: str, x_axis: tp.List[float], d_oss: tp.Dict[str, tp.Dict[str, float]], nome_outf: str, fpath: str) -> tp.NoReturn:
    '''Function to save observable(s) results '''
    file_data = open(fpath + os.sep + nome_outf, 'w')
    print(f"#L={L}", file = file_data)
    print(f"#extfield={h:.8f}", file = file_data)

    #Tries to make columns of datas
    fmt = f'#{x_name}\t\t'
    for oss in d_oss:
        fmt += '#{s}\t\t#d{s}\t\t'.format(s = oss)
    print(fmt.strip(), file = file_data)

    #Print datas on file, spaced with tab and with 8 digit precision
    for i, x in enumerate(x_axis):
        fmt = f'{x:.8f}\t'
        for oss in d_oss:
            fmt += ("{a[valore]["+ str(i) +"]:.8f}\t{a[errore]["+ str(i) +"]:.8f}\t").format(a = d_oss[oss])
        print(fmt.strip(),file=file_data)
    file_data.close()



def read_data(out_file: str) -> ising_type.tpoutdata:
    '''Function to read data from the output file of observable(s) calculation(s)'''
    file_data = open(out_file,'r')
    
    #Reads L, extfield
    L = int(sml.gread(file_data).strip().split('=')[1])
    h = float(sml.gread(file_data).strip().split('=')[1])
    
    #Reads the line with beta, observables and uncertainties
    s = sml.gread(file_data).strip().split()
    
    #Exclude the first charachter in file: #
    unitx = s[0][1:]
    oss_list = [ word[1:].lstrip('d') for word in s[1:] ]
    x_axis = []
    
    #Create the dictionary of values and errors for every observable in the file
    d_oss = { oss : {'valore': [], 'errore': []} for oss in oss_list }
    line = sml.gread(file_data).strip()
    
    while line:
        s = line.split()
        x_axis.append(float(s[0]))
        type = 'valore'
        for name, val in zip(oss_list, s[1:]):
            d_oss[name][type].append(float(val))
            type = type == 'valore' and 'errore' or 'valore'
        line = file_data.readline().strip()
    
    #Creating the dictionary with all the information loaded
    datas = {'x_axis': x_axis, 'd_oss': d_oss, 'L': L, 'extfield': h, 'unitx': unitx}
    return datas



def salva_storia(L: int, beta: float, extfield: float, vec: tp.Dict[str,tp.Dict[str, float]], file_name: str) -> tp.NoReturn:
    '''Function to save a single MC story '''

    file_data = open(file_name, 'w')
    print('#L=%d' %L, file=file_data)
    print('#beta=%.4f' %beta, file=file_data)
    print('#extfield=%.3f' %extfield, file=file_data)
    
    #Saving energy and/or magnetization
    for i in vec:
        if vec[i]:
            file_data.write( '#' + i + '=')
            for riga in vec[i]:
                file_data.write(str(riga)+ ' ')
            file_data.write('\n')
    file_data.close()






    
   
