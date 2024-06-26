#  A script to set the input parameters and run a CPMC calculation
#
# Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
# Package homepage: http://cpmc-lab.wm.edu
# Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
# Any publications resulting from either applying or building on the present package 
#   should cite the following journal article (in addition to the relevant literature on the method):

import numpy as np
import matplotlib.pyplot as plt
import CPMC_Lab as CP
from datetime import datetime 

## system parameters:
Lx=2; # The number of lattice sites in the x direction
Ly=1; # The number of lattice sites in the y direction
Lz=1; # The number of lattice sites in the z direction

N_up=1; # The number of spin-up electrons
N_dn=1; # The number of spin-down electrons

kx=+0.0819; # The x component of the twist angle in TABC (twist-averaging boundary condition)
ky=0; # The y component of the twist angle in TABC
kz=0; # The z component of the twist angle in TABC

U=4.0; # The on-site repulsion strength in the Hubbard Hamiltonian
tx=1; # The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; # The hopping amplitude between nearest neighbor sites in the y direction
tz=1; # The hopping amplitude between nearest neighbor sites in the z direction

## run parameters:
deltau=0.01; # The imaginary time step
N_wlk=100; # The number of random walkers
N_blksteps=40; # The number of random walk steps in each block
N_eqblk=2; # The number of blocks used to equilibrate the random walk before energy measurement takes place
N_blk=20; # The number of blocks used in the measurement phase
itv_modsvd=5; # The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers. No re-orthonormalization if itv_modsvd > N_blksteps
itv_pc=10; # The interval between two adjacent population controls. No population control if itv_pc > N_blksteps
itv_Em=20; # The interval between two adjacent energy measurements
now = datetime.now()
suffix= now.strftime("%d-%m-%Y-%H-%M-%S"); # time stamp for the saved *.mat filename. Can be changed to any desired string 

## invoke the main function
[E_ave,E_err,savedFile]=CP.CPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix);

print("E_ave",E_ave, "E_err",E_err)
## post-run:
# load saved data into workspace for post-run analysis:
data=np.load(savedFile+".npz");
# plot energy vs imaginary time
plt.figure();
plt.plot (N_blksteps*np.arange(0,N_blk,1)*deltau,data['E']);
#plt.axhline(y=E_ave,'g:')
plt.xlabel ('tau');
plt.ylabel ('E');
plt.show()
## Explanation of saved quantities:
# E: the array of energy of each block
# time: The total computational time
# E_nonint_v: the non-interacting energy levels of the system
# Phi_T: the trial wave function
# For other saved quantities, type "help CPMC_Lab"
