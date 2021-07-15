
import numpy as np
from numpy import linalg as LA
#import matplotlib
#import pylab
from scipy.sparse.linalg import eigsh
#import os


# define parameters

Iter =26

#h
hpoints = 50
hval = 0.0
hcenter=0.0
hin=-0.001
hfin=0.001


#m range
mval = 10
min=26
mfin=27

# parameters for my functions
gpoints = 111
gin = 0.4
gfin = 1.5
gcenter = 1.0

#to calculate magnetization
probe=10

# what if you don't want linearly spaced points for g? here I am
a = -(gcenter - gin) ** (1. / 3.)
b = (gfin - gcenter) ** (1. / 3.)
print(a, b)
gaxes = np.linspace(a, b, gpoints)

gcubic = 0.7 * gaxes ** 3 + gcenter * np.ones(gpoints) + 0.3 * (gaxes)
glinear = np.linspace(gin, gfin, gpoints)

#and also for h
ah = -(hcenter-hin) ** (1. / 3.)
bh = (hfin-hcenter) ** (1. / 3.)
print(ah, bh)
haxes = np.linspace(ah, bh, hpoints)

hcubic = 0.7 * haxes ** 3 + hcenter * np.ones(hpoints) + 0.3 * (haxes)
hlinear = np.linspace(hin, hfin, hpoints)
print(hcubic)
# define operatorsR
sigmax = np.array([[0, 1], [1, 0]])
sigmaz = np.array([[1, 0], [0, -1]])
ss = 'm' + str(min)+str(mfin) + 'iter' + str(Iter) + 'h' + str(hval)
# preparing to write on file
f = open("res2.txt", "w")  # file for the results
e = open("magnemma.txt", "w")  # for the magnetizations only


# this function gets the eigenvectors(psi) ordered by lowest eigenvalue(energy).
# how to return local variables?
def sorting(energy, psi):
    idx = np.argsort(energy)
    energy = energy[idx]
    psi = psi[:, idx]

for m in [mval]:
    for h in [hval]:

        for g in glinear:


            # define Block operators as the ones on single site
            Hsite = (-g) * sigmax - h * sigmaz
            BlockH = (-g) * sigmax - h * sigmaz
            Blocksigmaz = sigmaz
            Blocksigmax = sigmax
            I = np.eye(2)

            # DMRG
            for i in range(1, Iter):
                # control the progress
                print(m, h, g, i)

                # Building the Enlarged Block operators
                # He=np.kron(BlockH, np.eye(2))+np.kron(I,Hsite)-j*np.kron(Blocksigmaz, sigmaz)
                He = np.kron(BlockH, np.eye(2)) - g * np.kron(I, sigmax) - h * np.kron(I, sigmaz) - np.kron(Blocksigmaz,
                                                                                                            sigmaz)
                He = 0.5 * (He + He.T)
                Blocksigmaz = np.kron(I, sigmaz)
                Blocksigmax = np.kron(I, sigmax)
                I = np.kron(I, np.eye(2))

                if i > probe:
                    Blocksigmaz10 = np.kron(Blocksigmaz10, np.eye(2))
                    Blocksigmax10 = np.kron(Blocksigmax10, np.eye(2))

                # Building superblock operators
                Hsuper = np.kron(He, np.eye(len(He))) + np.kron(np.eye(len(He)), He) - np.kron(Blocksigmaz, Blocksigmaz)
                Hsuper = 0.5 * (Hsuper + Hsuper.T)

                # diagonalizing Hsuper and sorting eigenvectors. The sorting function was tested
                #energy, psi = LA.eigh(Hsuper)
                # sorting(energy,psi) #it does not return local variables for now, so it can be omitted
                energy, psi=eigsh(Hsuper,k=1,which='SA')
                idx = np.argsort(energy)
                energy = energy[idx]
                psi = psi[:, idx]

                # creating the density matrix for the superblock and taking the partial trace
                psi0 = psi[:, 0]
                nr = len(psi)
                Dim = int(np.sqrt(nr))
                PsiMatrix = psi0.reshape(Dim, Dim)
                Rho = np.dot(PsiMatrix, PsiMatrix)
                Rho = 0.5 * (Rho + Rho.T)

                # diagonalizing Rho and sorting eigenvectors
                val, vec = LA.eigh(Rho)
                idx = len(val) - 1 - np.argsort(val)
                val = val[idx]
                vec = vec[:, idx]

                # calculating the Von Neumann entropy of the block. build a function!
                entropy = 0
                for j in range(0, len(val)):  # riscrivi
                    if val[j] > 0:
                        entropy += -val[j] * np.log(val[j])

                # deciding how many states to keep and constructing the truncation operator
                #Nkeep = min(len(Rho), int(m))
                if(len(Rho)<m):
                    Nkeep=len(Rho)
                else:
                    Nkeep=m
                truncationdelta = 1 - np.sum(val[0:Nkeep])
                omatr = vec[:, 0:Nkeep]

                # calculating observables. Build a function!
                # if (i+1)<Iter:
                #if (i > 12) and ((i%10)==0):
                if (i==25 or i==35 or i==45 or i==60or i==80 or i==100):
                #if (i == 100 or i == 200 or i == 500 or i == 1000):
                #if(i==(Iter-1)):
                    a = psi0.reshape(len(psi0), 1)
                    b = psi0.reshape(1, len(psi0))
                    rhosuper = a @ b
                    sigmaz10super = np.kron(Blocksigmaz10, I)
                    sigmazL10super = np.kron(I, Blocksigmaz10)
                    sigmax10super = np.kron(Blocksigmax10, I)
                    sigmaxL10super = np.kron(I, Blocksigmax10)

                    Mz = np.sqrt(np.trace(sigmaz10super @ sigmazL10super @ rhosuper))
                    Mx = np.sqrt(np.trace(sigmax10super @ sigmaxL10super @ rhosuper))
                    # Mz=np.sqrt(np.trace(Blocksigmaz10 @ Blocksigmaz10 @ rhosuper))

                    d = 1
                    e.write("{}\t {} \t {} \t {} \t {} \t{} \t{} \t{}\n".format(d, i, g, h, m, Iter, Mx, Mz))

                # truncating the operators
                BlockH = omatr.T @ He @ omatr
                Blocksigmaz = omatr.T @ Blocksigmaz @ omatr
                Blocksigmax = omatr.T @ Blocksigmax @ omatr
                I = np.eye(Nkeep)

                if i == probe:  # to calculate correlations between distant sites
                    Blocksigmaz10 = Blocksigmaz
                    Blocksigmax10 = Blocksigmax

                if i > probe:
                    Blocksigmaz10 = omatr.T @ Blocksigmaz10 @ omatr
                    Blocksigmax10 = omatr.T @ Blocksigmax10 @ omatr

                # writing results on the file
                # if ((i%10==0)&(i>25)):
                if ((i % 300)== 0):
                    d = 1
                    f.write("{} \t{}\t {} \t{} \t {}\t {} \t {} \t {} \n".format(d, i, g, h, m, Iter, entropy, truncationdelta))

    # closing the file
f.close()
e.close()

print(ss)

# COMANDI UTIL


# E = [i for i in range(10) if i%2==0]
# print(E)

# np.dot, np.ein_sum, np.kron, np.linalg.eig
# dir(oggetto) help(oggetto)
#np.arange()
# ax, fig = plt.subplots() plt.grid() plt.legend() plt.plot(x,y, label="nome 1 ") #scatter  plt.show()
# ax.set_xlim ax.set_ylim

