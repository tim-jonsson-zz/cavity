import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

def readN():
    f1 = open("pulse.txt",'r')
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
    exNt = readN()
    fig,ax = plt.subplots(1,1)
    ax.plot(exNt[:,0],exNt[:,1])  
    ax.set_xlabel('t')
    plt.savefig(f"plots/pulse{name}.png")
    plt.clf()
    
#plotNt('2.0_S0a')
#plotNt('2.0_S0b')    
plotNt('_superG_c80_1')    
