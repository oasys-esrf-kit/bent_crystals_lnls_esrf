

import os
import numpy as np
import matplotlib.pyplot as plt

path = os.getcwd()+'\\'
fsize = 14

# Flux
filename1='res_APS_BLM_Si311_nrm5-1-1_sag01-1_ML_energy.txt'
filename2='res_APS_BLM_Si311_nrm5-1-1_sag01-1_PP_energy.txt'
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y = data1[:,2]
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y = data2[:,2]
data3_x = [30.,40.,50.,60.]
data3_y = [2.5e11,1.5e11,2.0e11,1.5e11]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Flux')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='PP')
ax1.plot(data3_x, data3_y, linestyle='-', color='green', linewidth=2, marker='o', label='paper')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('Flux [Ph/s/100mA]', fontsize=12)
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=2)
ax1.set_xscale('linear')
ax1.set_yscale('log')
plt.minorticks_on()
plt.tick_params(which='both', axis='both', direction='in', right=True, top=True, labelsize=12)
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'APS_BLM_Flux.png')
plt.close(fig)

# Flux density
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y = data1[:,3]
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y = data2[:,3]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Flux density')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='PP')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('Flux density [Ph/s/um²/100mA]', fontsize=12)
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=2)
ax1.set_xscale('linear')
ax1.set_yscale('log')
plt.minorticks_on()
plt.tick_params(which='both', axis='both', direction='in', right=True, top=True, labelsize=12)
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'APS_BLM_Fluxdensity.png')
plt.close(fig)

# beam size
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y1 = data1[:,4]
data1_y2 = data1[:,5]
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y1 = data2[:,4]
data2_y2 = data2[:,5]
data3_x2 = [30.,40.,50.,60.]
data3_y2 = [0.7,0.91,1.2,1.5]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Beam size')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='x ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='x PP')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='z ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='z PP')
ax2.plot(data3_x2, data3_y2, linestyle=':', color='green', linewidth=2, marker='o', label='z paper')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('BS X [mm]', fontsize=12)
ax2.set_ylabel('BS Z [mm]', fontsize=12)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
all_lines = lines1 + lines2
all_labels = labels1 + labels2
ax1.legend(all_lines, all_labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=3)
ax1.set_xscale('linear')
ax1.set_yscale('linear')
ax2.set_yscale('linear')
ax1.set_ylim([0.0, 0.5])
ax2.set_ylim([0.0, 3.0])
plt.minorticks_on()
for ax in (ax1, ax2):
    ax.tick_params(
        which='both',
        labelsize=12,
        direction='in'
    )
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'APS_BLM_BS.png')
plt.close(fig)

# beam divergence
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y1 = data1[:,6]
data1_y2 = data1[:,7]
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y1 = data2[:,6]
data2_y2 = data2[:,7]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Beam divergence')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='dx ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='dx PP')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='dz ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='dz PP')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('DIV X [mrad]', fontsize=12)
ax2.set_ylabel('DIV Z [mrad]', fontsize=12)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
all_lines = lines1 + lines2
all_labels = labels1 + labels2
ax1.legend(all_lines, all_labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=2)
ax1.set_xscale('linear')
ax1.set_yscale('linear')
ax2.set_yscale('linear')
ax1.set_ylim([1.0, 2.0])
ax2.set_ylim([0.0, 0.1])
plt.minorticks_on()
for ax in (ax1, ax2):
    ax.tick_params(
        which='both',
        labelsize=12,
        direction='in'
    )
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'APS_BLM_DIV.png')
plt.close(fig)

# energy resolution
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y = data1[:,8]/data1_x
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y = data2[:,8]/ data2_x
data3_x = [30.,40.,50.,60.]
data3_y = [9e-5,1.5e-4,2.4e-4,4.0e-4]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Energy resolution')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='PP')
ax1.plot(data3_x, data3_y, linestyle='-', color='green', linewidth=2, marker='o', label='paper')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('Resolution [10E-3]', fontsize=12)
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=3)
ax1.set_xscale('linear')
ax1.set_yscale('log')
plt.minorticks_on()
plt.tick_params(which='both', axis='both', direction='in', right=True, top=True, labelsize=12)
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'APS_BLM_resolution.png')
plt.close(fig)

# Laue Rs and Rm
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y1 = data1[:,10]*1e-3
data1_y2 = data1[:,11]*1e-3
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y1 = data2[:,10]*1e-3
data2_y2 = data2[:,11]*1e-3
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Rs Rm - Laue')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='Rs ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='Rs PP')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='PP')
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
plt.savefig(path+'APS_BLM_Rsm.png')
plt.close(fig)

