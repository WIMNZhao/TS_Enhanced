
import numpy as np
import math,sys
import pandas as pd
from rdkit import Chem

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


nprod = 376 * 500 * 500
top100 = {}
df = pd.read_csv("quinazoline_exhaustive.csv")

smiles = df['SMILES'].to_list()
for smi in smiles:
    mol = Chem.MolFromSmiles(smi)
    smi = Chem.MolToSmiles(mol)
    top100[smi] = None
#######################################################
df = pd.read_csv("test.csv")
smiles = df['SMILES'].to_list()
eff = df['Name'].to_list()

percent_lib = [(i+1)/nprod*100 for i in range(len(eff))]
recovery = []
sampling_eff = []
count = 0
for i in range(len(eff)):
    sampling_eff.append(float(eff[i].split('_')[-1]))
    if smiles[i] in top100:
       count += 1
    recovery.append(count)

print(count)

#==================================================================
fig, ax = plt.subplots(2,figsize=(6,8))
fig.align_labels()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.16)
plt.rcParams['axes.titlepad'] = 0 #-14  # pad is in points...
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

for i in range(2):
    #ax[i].spines['right'].set_visible(False)
    #ax[i].spines['top'].set_visible(False)

    ax[i].tick_params(size=6,width=1)
    ax[i].tick_params(which='minor',width=1, size=4)
    ax[i].tick_params(labelsize=12)
#==================================================================

ax[0].set_title('A',fontweight='bold', x=-0.13,y=0.94,fontsize=16) 
ax[0].plot(percent_lib,recovery)
ax[0].yaxis.set_major_locator(MultipleLocator(20))
ax[0].yaxis.set_minor_locator(MultipleLocator(10))
#ax[0].xaxis.set_major_locator(MultipleLocator(0.05))
#ax[0].xaxis.set_minor_locator(MultipleLocator(0.025))
ax[0].set_ylim(-10,110)
ax[0].set_ylabel(r'recovery (%)', fontsize=14)

ax[1].set_title('B',fontweight='bold', x=-0.13,y=0.94,fontsize=16) 
ax[1].plot(percent_lib,sampling_eff)
ax[1].yaxis.set_major_locator(MultipleLocator(20))
ax[1].yaxis.set_minor_locator(MultipleLocator(10))
#ax[1].xaxis.set_major_locator(MultipleLocator(0.05))
#ax[1].xaxis.set_minor_locator(MultipleLocator(0.025))
ax[1].set_xlabel(r'percent of library screened (%)', fontsize=14)
ax[1].set_ylabel(r'sampling efficiency (%)', fontsize=14)
ax[1].set_ylim(10,110)

#==========================================================

plt.savefig('test.png',bbox_inches='tight',dpi=300) #,transparent=True)
