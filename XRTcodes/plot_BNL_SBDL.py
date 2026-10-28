

import os
import numpy as np
import matplotlib.pyplot as plt

path = os.getcwd()+'\\'
fsize = 14

# Flux
filename1='res_BNL_SBDL_highfluxSi111_nrm100_sag0-11_ML_energy.txt'
filename2='res_BNL_SBDL_highfluxSi111_nrm100_sag0-11_PP_energy.txt'
filename3='res_BNL_SBDL_highresSi111_nrm100_sag0-11_ML_energy.txt'
filename4='res_BNL_SBDL_highresSi111_nrm100_sag0-11_PP_energy.txt'
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
data5_y = [3.0e12,4.75e12,4.25e12,2.23e12,1.1e12]
data6_y = [6.5e11,9.25e11,4.5e11,1.73e11,5.0e10]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Flux')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='hFlux-ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='hFlux-PP')
ax1.plot(data3_x, data3_y, linestyle='-', color='red', linewidth=2, marker='o', label='hRes-ML')
ax1.plot(data4_x, data4_y, linestyle='-', color='red', linewidth=2, marker='s', label='hRes-PP')
ax1.plot(data1_x, data5_y, linestyle='-', color='green', linewidth=2, marker='o', label='paper-hFlux')
ax1.plot(data1_x, data6_y, linestyle='-', color='green', linewidth=2, marker='s', label='paper-hRes')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('Flux [Ph/s/100mA]', fontsize=12)
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=3)
ax1.set_xscale('linear')
ax1.set_yscale('log')
plt.minorticks_on()
plt.tick_params(which='both', axis='both', direction='in', right=True, top=True, labelsize=12)
plt.grid(which='both', alpha=0.2)    
plt.tight_layout()
plt.savefig(path+'BNL_SBDL_Flux.png')
plt.close(fig)

# Flux density
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y = data1[:,3]
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y = data2[:,3]
data3 = np.genfromtxt(path+filename3, skip_header=1)
data3_x = data3[:,0]*1e-3
data3_y = data3[:,3]
data4 = np.genfromtxt(path+filename4, skip_header=1)
data4_x = data4[:,0]*1e-3
data4_y = data4[:,3]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Flux density')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='hFlux-ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='hFlux-PP')
ax1.plot(data3_x, data3_y, linestyle='-', color='red', linewidth=2, marker='o', label='hFlux-ML')
ax1.plot(data4_x, data4_y, linestyle='-', color='red', linewidth=2, marker='s', label='hFlux-PP')
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
plt.savefig(path+'BNL_SBDL_Fluxdensity.png')
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
data3 = np.genfromtxt(path+filename3, skip_header=1)
data3_x = data3[:,0]*1e-3
data3_y1 = data3[:,4]
data3_y2 = data3[:,5]
data4 = np.genfromtxt(path+filename4, skip_header=1)
data4_x = data4[:,0]*1e-3
data4_y1 = data4[:,4]
data4_y2 = data4[:,5]
data5_y1 = [0.71,0.81,0.61,0.64,0.78]
data6_y1 = [0.08,0.06,0.057,0.055,0.064]
data5_y2 = [0.65,0.65,0.80,0.72,0.64]
data6_y2 = [2.6,2.3,2.0,1.7,1.5]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Beam size')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='x hFlux-ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='x hFlux-PP')
ax1.plot(data3_x, data3_y1, linestyle='-', color='red', linewidth=2, marker='o', label='x hRes-ML')
ax1.plot(data4_x, data4_y1, linestyle='-', color='red', linewidth=2, marker='s', label='x hRes-PP')
ax1.plot(data1_x, data5_y1, linestyle='-', color='green', linewidth=2, marker='o', label='paper x-hFlux')
ax1.plot(data1_x, data5_y2, linestyle='-', color='green', linewidth=2, marker='s', label='paper x-hRes')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='z hFlux-ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='z hFlux-PP')
ax2.plot(data3_x, data3_y2, linestyle=':', color='red', linewidth=2, marker='o',label='z hRes-ML')
ax2.plot(data4_x, data4_y2, linestyle=':', color='red', linewidth=2, marker='s',label='z hRes-PP')
ax2.plot(data1_x, data6_y1, linestyle=':', color='green', linewidth=2, marker='o', label='paper z-hFlux')
ax2.plot(data1_x, data6_y2, linestyle=':', color='green', linewidth=2, marker='s', label='paper z-hRes')
ax1.set_xlabel('Energy [keV]', fontsize=12)
ax1.set_ylabel('BS X [mm]', fontsize=12)
ax2.set_ylabel('BS Z [mm]', fontsize=12)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
all_lines = lines1 + lines2
all_labels = labels1 + labels2
ax1.legend(all_lines, all_labels, bbox_to_anchor=(0.5, -0.15), loc="upper center", fancybox=True, shadow=True, ncol=6)
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
plt.savefig(path+'BNL_SBDL_BS.png')
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
data3 = np.genfromtxt(path+filename3, skip_header=1)
data3_x = data3[:,0]*1e-3
data3_y1 = data3[:,6]
data3_y2 = data3[:,7]
data4 = np.genfromtxt(path+filename4, skip_header=1)
data4_x = data4[:,0]*1e-3
data4_y1 = data4[:,6]
data4_y2 = data4[:,7]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Beam divergence')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='dx hFlux-ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='dx hFlux-PP')
ax1.plot(data3_x, data3_y1, linestyle='-', color='red', linewidth=2, marker='o', label='dx hRes-ML')
ax1.plot(data4_x, data4_y1, linestyle='-', color='red', linewidth=2, marker='s', label='dx hRes-PP')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='dz hFlux-ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='dz hFlux-PP')
ax2.plot(data3_x, data3_y2, linestyle=':', color='red', linewidth=2, marker='o',label='dz hRes-ML')
ax2.plot(data4_x, data4_y2, linestyle=':', color='red', linewidth=2, marker='s',label='dz hRes-PP')
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
plt.savefig(path+'BNL_SBDL_DIV.png')
plt.close(fig)

# energy resolution
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y = data1[:,8]/data1_x
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y = data2[:,8]/ data2_x
data3 = np.genfromtxt(path+filename3, skip_header=1)
data3_x = data3[:,0]*1e-3
data3_y = data3[:,8]/data3_x
data4 = np.genfromtxt(path+filename4, skip_header=1)
data4_x = data4[:,0]*1e-3
data4_y = data4[:,8]/data4_x
data5_y = [4.0e-4,5.9e-4,9.9e-4,1.6e-3,2.3e-3]
data6_y = [1.0e-4,1.2e-4,1.2e-4,1.3e-4,1.4e-4]
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Energy resolution')
ax1.plot(data1_x, data1_y, linestyle='-', color='blue', linewidth=2, marker='o', label='hFlux-ML')
ax1.plot(data2_x, data2_y, linestyle='-', color='blue', linewidth=2, marker='s', label='hFlux-PP')
ax1.plot(data3_x, data3_y, linestyle='-', color='red', linewidth=2, marker='o', label='hRes-ML')
ax1.plot(data4_x, data4_y, linestyle='-', color='red', linewidth=2, marker='s', label='hRes-PP')
ax1.plot(data1_x, data5_y, linestyle='-', color='green', linewidth=2, marker='o', label='paper-hFlux')
ax1.plot(data1_x, data6_y, linestyle='-', color='green', linewidth=2, marker='s', label='paper-hRes')
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
plt.savefig(path+'BNL_SBDL_resolution.png')
plt.close(fig)

# Laue Rs and Rm
data1 = np.genfromtxt(path+filename1, skip_header=1)
data1_x = data1[:,0]*1e-3
data1_y1 = data1[:,9]*1e-3
data1_y2 = data1[:,10]*1e-3
data2 = np.genfromtxt(path+filename2, skip_header=1)
data2_x = data2[:,0]*1e-3
data2_y1 = data2[:,9]*1e-3
data2_y2 = data2[:,10]*1e-3
data3 = np.genfromtxt(path+filename3, skip_header=1)
data3_x = data3[:,0]*1e-3
data3_y1 = data3[:,9]*1e-3
data3_y2 = data3[:,10]*1e-3
data4 = np.genfromtxt(path+filename4, skip_header=1)
data4_x = data4[:,0]*1e-3
data4_y1 = data4[:,9]*1e-3
data4_y2 = data4[:,10]*1e-3
fig, ax1 = plt.subplots(figsize=(10,6))
plt.title('Rs Rm - Laue')
ax1.plot(data1_x, data1_y1, linestyle='-', color='blue', linewidth=2, marker='o', label='Rs hFlux-ML')
ax1.plot(data2_x, data2_y1, linestyle='-', color='blue', linewidth=2, marker='s', label='Rs hFlux-PP')
ax1.plot(data3_x, data3_y1, linestyle='-', color='red', linewidth=2, marker='o', label='Rs hRes-ML')
ax1.plot(data4_x, data4_y1, linestyle='-', color='red', linewidth=2, marker='s', label='Rs hRes-PP')
ax2=plt.twinx()
ax2.plot(data1_x, data1_y2, linestyle=':', color='blue', linewidth=2, marker='o', label='Rm hFlux-ML')
ax2.plot(data2_x, data2_y2, linestyle=':', color='blue', linewidth=2, marker='s',label='Rm hFlux-PP')
ax2.plot(data3_x, data3_y2, linestyle=':', color='red', linewidth=2, marker='o',label='Rm hRes-ML')
ax2.plot(data4_x, data4_y2, linestyle=':', color='red', linewidth=2, marker='s',label='Rm hRes-PP')
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
plt.savefig(path+'BNL_SBDL_Rsm.png')
plt.close(fig)

