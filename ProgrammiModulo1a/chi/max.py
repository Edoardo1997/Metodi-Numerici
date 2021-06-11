
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pylab
from scipy.optimize import curve_fit
def ff(x, aa, bb,cc):
    return aa-(bb*(x-cc)**2)
L_list=list(range(10,52,5))
data={L: pylab.loadtxt('chiL'+str(L)+'.txt',unpack = True) for L in L_list}
estr={10: (0.38,0.42), 15: (0.40,0.42), 20: (0.413,0.434), 25: (0.415,0.438), 30: (0.4185,0.439), 35: (0.424,0.4345), 40: (0.424,0.439), 45: (0.428,0.435), 50: (0.428,0.4355), 60: (0.427,0.4365)}
datan={L: [[],[],[]] for L in L_list}
for L in L_list:
    for i in range(len(data[L][0])):
        if estr[L][0]<=data[L][0][i]<=estr[L][1]:
            datan[L][0].append(data[L][0][i])
            datan[L][1].append(data[L][1][i])
            datan[L][2].append(data[L][2][i])
    

f = open("suscettivitamax.txt", "w")
print('L    beta_pc     Dbeta_pc    chi_max  Dchi_max   Chi/ndf',file=f)
for L in L_list:
    init=(L,100, 0.43)
    sigma=np.array(datan[L][2])
    x=np.array(datan[L][0])
    y=np.array(datan[L][1])
    w=1/sigma**2
    #fit
    pars,covm=curve_fit(ff,x,y,init,sigma, absolute_sigma=False)

    chi2 = ((w*(y-ff(x,*pars))**2)).sum()


    # determine the ndof
    ndof=len(x)-2
    print('Lunghezza del reticolo = ', L)
    print ('chi2=',chi2, 'ndof=',ndof)
    print('Il punto di max è', pars[2])
    print('Il max è', pars[0],'\n')
    errorepuntomax = np.sqrt(covm[2,2])
    erroremax = np.sqrt(covm[0,0])
    cr=chi2/ndof

    datifile = str(L) + '   ' + str(pars[2]) + '   '+  str(errorepuntomax) + '   '+ str(pars[0])+ '   '+ str(erroremax)+'  '+str(cr)+'\n'

    f.write(datifile)

f.close()
'''
x0, y0, dy0 = pylab.loadtxt('chiL10.txt',unpack = True)
x1, y1, dy1 = pylab.loadtxt('chiL20.txt',unpack = True)
x2, y2, dy2 = pylab.loadtxt('chiL30.txt',unpack = True)
x3, y3, dy3 = pylab.loadtxt('chiL40.txt',unpack = True)
x4, y4, dy4 = pylab.loadtxt('chiL50.txt',unpack = True)
x5, y5, dy5 = pylab.loadtxt('chiL60.txt',unpack = True)

x0n=[]
y0n=[]
d0n=[]
for i in x0:
    if 0.40 <=i<0.43:
        x0n.append(i)
        y0n.append(i)
        y0n.append(i)

# Creo il file in cui mettere i dati L, max della suscettività, punto di max


def fit(L,x,y,dy,first,end):

    file_one = open("suscettivitamax.txt", "a")

    x=x[first:end]
    y=y[first:end]
    dy=dy[first:end]

   # pylab.errorbar(x,y,dy,linestyle = '', color = 'black', marker = '.')



    #fit
    sigma=dy
    w=1/sigma**2

    def ff(x, aa, bb,cc):
        return aa-(bb*(x-cc)**2)

    init=(L,100, 0.43)

    pars,covm=curve_fit(ff,x,y,init,sigma, absolute_sigma=False)

    chi2 = ((w*(y-ff(x,*pars))**2)).sum()


    # determine the ndof
    ndof=len(x)-2
    print('Lunghezza del reticolo = ', L)
    print ('chi2=',chi2, 'ndof=',ndof)
    print('Il punto di max è', pars[2])
    print('Il max è', pars[0],'\n')



    errorepuntomax = np.sqrt(covm[2,2])
    erroremax = np.sqrt(covm[0,0])


    datifile = str(L) + '   ' + str(pars[2]) + '   '+  str(errorepuntomax) + '   '+ str(pars[0])+ '   '+ str(erroremax)

    file_one.write(datifile)
    file_one.write("\n")
    file_one.close()

    return(pars,covm)




L0=10
L1=20
L2=30
L3=40
L4=50
L5=60

file_one = open("suscettivitamax.txt", "w")
print('L beta_pc Dbeta_pc chi_max Dchi_max',file=file_one)
file_one.close()
''''''
m0,s0=fit(L0,x0,y0,dy0,105,140)
m1,s1=fit(L1,x1,y1,dy1,60,170)
m2,s2=fit(L2,x2,y2,dy2,90,165)
m3,s3=fit(L3,x3,y3,dy3,100,138)
m4,s4=fit(L4,x4,y4,dy4,105,140)
m5,s5=fit(L5,x5,y5,dy5,105,140)


s1[0][0]=np.sqrt(s1[0][0])
s2[0][0]=np.sqrt(s2[0][0])
s3[0][0]=np.sqrt(s3[0][0])
s4[0][0]=np.sqrt(s4[0][0])

gam= -7/4 #gamma su nu

Y=(m1[0]*L1**gam,m2[0]*L2**gam,m3[0]*L3**gam,m4[0]*L4**gam)
X=(L1,L2,L3,L4)
DY=(s1[0][0]*L1**gam,s2[0][0]*L2**gam,s3[0][0]*L3**gam,s4[0][0]*L4**gam)

#print(np.log(m1[0]/m2[0])/np.log(20/30))
pylab.errorbar(X,Y,DY,linestyle = '', color = 'black', marker = 'o')

'''







