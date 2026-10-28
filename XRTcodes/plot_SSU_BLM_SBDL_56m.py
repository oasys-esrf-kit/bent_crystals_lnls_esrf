

import os
import numpy as np
import matplotlib.pyplot as plt

path = os.getcwd()+'\\'
fsize = 14

# Flux
filename1='res_SSU_BLM_56m_Si311_nrm1911_sag01-1_ML_energy.txt'
filename2='res_SSU_BLM_56m_Si311_nrm1911_sag01-1_PP_energy.txt'
filename3='res_SSU_SBDL_56m_Si111_nrm100_sag011_ML_energy.txt'
filename4='res_SSU_SBDL_56m_Si111_nrm100_sag011_PP_energy.txt'
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y = data1[:,2]
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y = data2[:,2]
data3 = np.genfromtxt(path+filename3, skip_header=1)
data3_x = data3[:,0]*1e-3
data3_y = data3[:,2]
data4 = np.genfromtxt(path+filename4, skip_header=1)
data4_x = data4[:,0]*1e-3
data4_y = data4[:,2]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Flux')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='BLM-ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='BLM-PP')
ax1.plot(data3_x, data3_y, linestyle='-', color='red', linewidth=2, marker='o', label='SBDL-ML')
ax1.plot(data4_x, data4_y, linestyle='-', color='red', linewidth=2, marker='s', label='SBDL-PP')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('Flux [Ph/s/100mA]', fontsize=12)
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=4)
ax1.set_xscale('linear')
ax1.set_yscale('log')
plt.minorticks_on()
plt.tick_params(which='both', axis='both', direction='in', right=True, top=True, labelsize=12)
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'SSU_56m_Flux.png')
plt.close(fig)

# Flux density
data1_y = data1[:,3]
data2_y = data2[:,3]
data3_y = data3[:,3]
data4_y = data4[:,3]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Flux density')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='BLM-ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='BLM-PP')
ax1.plot(data3_x, data3_y, linestyle='-', color='red', linewidth=2, marker='o', label='SBDL-ML')
ax1.plot(data4_x, data4_y, linestyle='-', color='red', linewidth=2, marker='s', label='SBDL-PP')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('Flux density [Ph/s/um²/100mA]', fontsize=12)
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=4)
ax1.set_xscale('linear')
ax1.set_yscale('log')
plt.minorticks_on()
plt.tick_params(which='both', axis='both', direction='in', right=True, top=True, labelsize=12)
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'SSU_56m_Fluxdensity.png')
plt.close(fig)

# beam size
data1_y1 = data1[:,4]
data1_y2 = data1[:,5]
data2_y1 = data2[:,4]
data2_y2 = data2[:,5]
data3_y1 = data3[:,4]
data3_y2 = data3[:,5]
data4_y1 = data4[:,4]
data4_y2 = data4[:,5]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Beam size')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='x BLM-ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='x BLM-PP')
ax1.plot(data3_x, data3_y1, linestyle='-', color='red', linewidth=2, marker='o', label='x SBDL-ML')
ax1.plot(data4_x, data4_y1, linestyle='-', color='red', linewidth=2, marker='s', label='x SBDL-PP')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='z BLM-ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='z BLM-PP')
ax2.plot(data3_x, data3_y2, linestyle=':', color='red', linewidth=2, marker='o',label='z SBDL-ML')
ax2.plot(data4_x, data4_y2, linestyle=':', color='red', linewidth=2, marker='s',label='z SBDL-PP')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('BS X [mm]', fontsize=12)
ax2.set_ylabel('BS Z [mm]', fontsize=12)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
all_lines = lines1 + lines2
all_labels = labels1 + labels2
ax1.legend(all_lines, all_labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=4)
ax1.set_xscale('linear')
ax1.set_yscale('linear')
ax2.set_yscale('linear')
plt.minorticks_on()
for ax in (ax1, ax2):
    ax.tick_params(
        which='both',
        labelsize=12,
        direction='in'
    )
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'SSU_56m_BS.png')
plt.close(fig)

# beam divergence
data1_y1 = data1[:,6]
data1_y2 = data1[:,7]
data2_y1 = data2[:,6]
data2_y2 = data2[:,7]
data3_y1 = data3[:,6]
data3_y2 = data3[:,7]
data4_y1 = data4[:,6]
data4_y2 = data4[:,7]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Beam divergence')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='dx BLM-ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='dx BLM-PP')
ax1.plot(data3_x, data3_y1, linestyle='-', color='red', linewidth=2, marker='o', label='dx SBDL-ML')
ax1.plot(data4_x, data4_y1, linestyle='-', color='red', linewidth=2, marker='s', label='dx SBDL-PP')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='dz BLM-ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='dz BLM-PP')
ax2.plot(data3_x, data3_y2, linestyle=':', color='red', linewidth=2, marker='o',label='dz SBDL-ML')
ax2.plot(data4_x, data4_y2, linestyle=':', color='red', linewidth=2, marker='s',label='dz SBDL-PP')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('DIV X [mrad]', fontsize=12)
ax2.set_ylabel('DIV Z [mrad]', fontsize=12)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
all_lines = lines1 + lines2
all_labels = labels1 + labels2
ax1.legend(all_lines, all_labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=4)
ax1.set_xscale('linear')
ax1.set_yscale('linear')
ax2.set_yscale('linear')
plt.minorticks_on()
for ax in (ax1, ax2):
    ax.tick_params(
        which='both',
        labelsize=12,
        direction='in'
    )
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'SSU_56m_DIV.png')
plt.close(fig)

# energy resolution
data1_y = data1[:,8]/data1_x
data2_y = data2[:,8]/ data2_x
data3_y = data3[:,8]/data3_x
data4_y = data4[:,8]/data4_x
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Energy resolution')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='BLM-ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='BLM-PP')
ax1.plot(data3_x, data3_y, linestyle='-', color='red', linewidth=2, marker='o', label='SBDL-ML')
ax1.plot(data4_x, data4_y, linestyle='-', color='red', linewidth=2, marker='s', label='SBDL-PP')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('Resolution [10E-3]', fontsize=12)
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=4)
ax1.set_xscale('linear')
ax1.set_yscale('log')
plt.minorticks_on()
plt.tick_params(which='both', axis='both', direction='in', right=True, top=True, labelsize=12)
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'SSU_56m_resolution.png')
plt.close(fig)

# Laue Rs and Rm
data1_y1 = data1[:,10]*1e-3
data1_y2 = data1[:,11]*1e-3
data2_y1 = data2[:,10]*1e-3
data2_y2 = data2[:,11]*1e-3
data3_y1 = data3[:,9]*1e-3
data3_y2 = data3[:,10]*1e-3
data4_y1 = data4[:,9]*1e-3
data4_y2 = data4[:,10]*1e-3
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Rs Rm - Laue')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='Rs BLM-ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='Rs BLM-PP')
ax1.plot(data3_x, data3_y1, linestyle='-', color='red', linewidth=2, marker='o', label='Rs SBDL-ML')
ax1.plot(data4_x, data4_y1, linestyle='-', color='red', linewidth=2, marker='s', label='Rs SBDL-PP')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='Rm BLM-ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='Rm BLM-PP')
ax2.plot(data3_x, data3_y2, linestyle=':', color='red', linewidth=2, marker='o',label='Rm SBDL-ML')
ax2.plot(data4_x, data4_y2, linestyle=':', color='red', linewidth=2, marker='s',label='Rm SBDL-PP')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('Rs [m]', fontsize=12)
ax2.set_ylabel('Rm [m]', fontsize=12)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
all_lines = lines1 + lines2
all_labels = labels1 + labels2
ax1.legend(all_lines, all_labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=4)
ax1.set_xscale('linear')
ax1.set_yscale('linear')
ax2.set_yscale('linear')
plt.minorticks_on()
for ax in (ax1, ax2):
    ax.tick_params(
        which='both',
        labelsize=12,
        direction='in'
    )
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'SSU_56m_Rsm.png')
plt.close(fig)

