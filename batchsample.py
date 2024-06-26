# A script to loop over multiple sets of input parameters and run a CPMC calculation for each set

import numpy as np
import matplotlib.pyplot as plt
import CPMC_Lab as CP  

## system parameters:
Lx=2; # The number of lattice sites in the x direction
Ly=1; # The number of lattice sites in the y direction
Lz=1; # The number of lattice sites in the z direction

N_up=1; # The number of spin-up electrons
N_dn=1; # The number of spin-down electrons

kx=0; # The x component of the twist angle in TABC (twist-averaging boundary condition)
ky=0; # The y component of the twist angle in TABC
kz=0; # The z component of the twist angle in TABC

U=4.0; # The on-site repulsion strength in the Hubbard Hamiltonian
tx=1; # The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; # The hopping amplitude between nearest neighbor sites in the y direction
tz=1; # The hopping amplitude between nearest neighbor sites in the z direction

## run parameters:
deltau=0.01; # The imaginary time step
N_wlk=np.array([100,200,500]); # The number of random walkers
N_blksteps=40; # The number of random walk steps in each block
N_eqblk=2; # The number of blocks used to equilibrate the random walk before energy measurement takes place
N_blk=np.array([10,20,30]); # The number of blocks used in the measurement phase
itv_modsvd=5; # The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers. No re-orthonormalization if itv_modsvd > N_blksteps
itv_pc=10; # The interval between two adjacent population controls. No population control if itv_pc > N_blksteps
itv_Em=20; # The interval between two adjacent energy measurements

## Initialize the batch run
N_run=N_wlk.shape[0]; #replace argument by the parameter that needs to be looped over
E_ave=np.zeros((N_run,1));
E_err=np.zeros((N_run,1));

## invoke the main function
for i in range(0,N_run):
    suffix='_Nwlk'+str(N_wlk[i]); # Set the suffix to distinguish between different runs in the same batch    
    # call main function AFQMC_Hub
    [E_ave[i],E_err[i],savedFile]=CP.CPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk[i],N_blksteps,N_eqblk,N_blk[i],itv_modsvd,itv_pc,itv_Em,suffix);    

## post-run:
# plot energy vs different run parameters
#figure;
plt.figure()
plt.errorbar(N_wlk,E_ave,yerr=E_err);
plt.xlabel ('$N_{wlk}$');
plt.ylabel ('Energy');

## Explanation of saved quantities:
# E: the array of energy of each block
# time: The total computational time
# E_nonint_v: the non-interacting energy levels of the system
# Phi_T: the trial wave function
# For other saved quantities, type "help CPMC_Lab"
