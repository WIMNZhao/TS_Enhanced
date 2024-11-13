
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
nr = 10
smiles_all = []
eff_all = []
for i in range(0, nr):
    df = pd.read_csv(f"ts_result_1_{i:03d}.csv")
    smiles_all.append(df['SMILES'].to_list())
    eff_all.append(df['Name'].to_list())

leng = [len(x) for x in smiles_all]
mins = min(leng)

percent_lib = [(i+1)/nprod*100 for i in range(mins)]
recovery_all = []
sampling_eff_all = []

for smiles, eff in zip(smiles_all, eff_all):
    count = 0
    recovery = []
    sampling_eff = []
    for i in range(mins):
        sampling_eff.append(float(eff[i].split('_')[-1]))
        if smiles[i] in top100:
           count += 1
        recovery.append(count)
    recovery_all.append(recovery)
    sampling_eff_all.append(sampling_eff)

recovery_average = []
sampling_eff_average = []
recovery_std = []
sampling_eff_std = []
recovery_all = np.array(recovery_all)
sampling_eff_all = np.array(sampling_eff_all)
for i in range(mins):
    recovery_average.append(np.average(recovery_all[:,i]))
    recovery_std.append(np.std(recovery_all[:,i]))
    sampling_eff_average.append(np.average(sampling_eff_all[:,i]))
    sampling_eff_std.append(np.std(sampling_eff_all[:,i]))
r_aver = np.array(recovery_average)
r_std = np.array(recovery_std)
s_aver = np.array(sampling_eff_average)
s_std = np.array(sampling_eff_std)
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

ax[0].set_title('A',fontweight='bold', x=-0.14,y=0.94,fontsize=16) 
ax[0].plot(percent_lib,recovery_average,label=r'n_per_search_cycle = 1')
ax[0].fill_between(percent_lib,r_aver-r_std,r_aver+r_std,alpha=0.5)
ax[0].yaxis.set_major_locator(MultipleLocator(20))
ax[0].yaxis.set_minor_locator(MultipleLocator(10))
ax[0].xaxis.set_major_locator(MultipleLocator(0.05))
ax[0].xaxis.set_minor_locator(MultipleLocator(0.025))
ax[0].set_ylim(-10,110)
ax[0].set_xlim(-0.01,0.2)
ax[0].set_ylabel(r'recovery (%)', fontsize=14)

ax[1].set_title('B',fontweight='bold', x=-0.14,y=0.94,fontsize=16) 
ax[1].plot(percent_lib,sampling_eff_average)
ax[1].fill_between(percent_lib,s_aver-s_std,s_aver+s_std,alpha=0.5)
ax[1].yaxis.set_major_locator(MultipleLocator(20))
ax[1].yaxis.set_minor_locator(MultipleLocator(10))
ax[1].xaxis.set_major_locator(MultipleLocator(0.05))
ax[1].xaxis.set_minor_locator(MultipleLocator(0.025))
ax[1].set_xlabel(r'percent of library screened (%)', fontsize=14)
ax[1].set_ylabel(r'sampling efficiency (%)', fontsize=14)
ax[1].set_ylim(10,110)
ax[1].set_xlim(-0.01,0.2)

#==========================================================
nr = 10
smiles_all = []
eff_all = []
for i in range(0, nr):
    df = pd.read_csv(f"ts_result_3_{i:03d}.csv")
    smiles_all.append(df['SMILES'].to_list())
    eff_all.append(df['Name'].to_list())

leng = [len(x) for x in smiles_all]
mins = min(leng)

percent_lib = [(i+1)/nprod*100 for i in range(mins)]
recovery_all = []
sampling_eff_all = []

for smiles, eff in zip(smiles_all, eff_all):
    count = 0
    recovery = []
    sampling_eff = []
    for i in range(mins):
        sampling_eff.append(float(eff[i].split('_')[-1]))
        if smiles[i] in top100:
           count += 1
        recovery.append(count)
    recovery_all.append(recovery)
    sampling_eff_all.append(sampling_eff)

recovery_average = []
sampling_eff_average = []
recovery_std = []
sampling_eff_std = []
recovery_all = np.array(recovery_all)
sampling_eff_all = np.array(sampling_eff_all)
for i in range(mins):
    recovery_average.append(np.average(recovery_all[:,i]))
    recovery_std.append(np.std(recovery_all[:,i]))
    sampling_eff_average.append(np.average(sampling_eff_all[:,i]))
    sampling_eff_std.append(np.std(sampling_eff_all[:,i]))
r_aver = np.array(recovery_average)
r_std = np.array(recovery_std)
s_aver = np.array(sampling_eff_average)
s_std = np.array(sampling_eff_std)

ax[0].plot(percent_lib,recovery_average,label=r'n_per_search_cycle = 3')
ax[0].fill_between(percent_lib,r_aver-r_std,r_aver+r_std,alpha=0.5)
ax[1].plot(percent_lib,sampling_eff_average)
ax[1].fill_between(percent_lib,s_aver-s_std,s_aver+s_std,alpha=0.5)
ax[0].legend(frameon=False,loc='best',fontsize=12)

plt.savefig('fig4.png',bbox_inches='tight',dpi=300) #,transparent=True)
