import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

def readN(fname):
    f1 = open(fname,'r')
    lines = f1.readlines()
    f1.close()
    N_t = []
    for line in lines:
        exN_s = line.split()
        exN = []
        for ni in exN_s:
            exN.append(float(ni))
        N_t.append(exN)
    return np.array(N_t)
  
def plotNt(name):
    #fname = f'data/w{name}.txt'
    exNt = readN("data/Nin_Nout_lanc.txt")
    Nt = 2500
    tt = np.arange(1,Nt+1)*0.05
    fig,ax = plt.subplots(3,1)
    ax[0].plot(tt,exNt[:,0],label='n1')  
    ax[0].plot(tt,exNt[:,1],label='n2')  
    ax[0].set_ylim([-0.05,1.05])
    #ax[0].set_ylim([0.3,0.7])
    ax[1].plot(tt,exNt[:,2],label='n3')  
    ax[1].plot(tt,exNt[:,3],label='n4')  
    ax[1].set_ylim([-0.05,1.05])
    #ax[1].set_ylim([0.3,0.7])
    ax[2].plot(tt,exNt[:,4],label='np')  
    for i in range(3):
        ax[i].legend()
    plt.savefig(f"plots/exNt_w{name}.png")
    plt.clf()
    
#plotNt('2.0_S0a')
#plotNt('2.0_S0b')    
plotNt('1.0_DaSwoff_gpsuperGsin1_2.0_1')    
