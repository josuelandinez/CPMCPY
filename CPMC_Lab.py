import numpy as np
import time as tm
import pop_cntrl as pc
import stblz as stb
import stepwlk as stpw
import initialization as init

def CPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix):
    # function [E_ave,E_err,savedFileName]=CPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em, suffix)
    # Perform a constrained path Monte Carlo calculatiion. Main function in the CPMC-Lab package
    # Input
    #   Lx: The number of lattice sites in the x direction.
    #   Ly: The number of lattice sites in the y direction.
    #   Lz: The number of lattice sites in the z direction.
    #   N_up: The number of spin-up electrons
    #   N_dn: The number of spin-down electrons
    #   kx: The x component of the twist angle in TABC (twist-averaging boundary condition)
    #   ky: The y component of the twist angle in TABC
    #   kz: The z component of the twist angle in TABC
    #   U: The on-site repulsion strength in the Hubbard Hamiltonian
    #   tx: The hopping amplitude between nearest-neighbor sites in the x direction
    #   ty: The hopping amplitude between nearest neighbor sites in the y direction
    #   tz: The hopping amplitude between nearest neighbor sites in the z direction
    #   deltau: The imaginary time step
    #   N_wlk: The number of random walkers
    #   N_blksteps: The number of random walk steps in each block
    #   N_eqblk: The number of blocks used to equilibrate the random walk before energy measurement takes place
    #   N_blk: The number of blocks used in the measurement phase
    #   itv_modsvd: The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers.
    #   itv_pc: The interval between two adjacent population controls
    #   itv_Em: The interval between two adjacent energy measurements
    #   suffix: an identifying string e.g. timestamp to be appended to the end of the saved *.mat file
    # Output:
    #   E_ave: the ground state energy
    #   E_err: the standard error in the ground state energy
    #   savedFileName: name of the saved data file
    #
    
    ## Initialization
    start_time = tm.time(); # start the  timer
    [N_sites,N_par,H_k,Proj_k_half,E_nonint_v,Phi_T,E_T, Phi,w,O,fac_norm,aux_fld,savedFileName]=init.initialization(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix); # initialize internal constants, form the trial wave function and assemble the initial population of walkers
    #format long;
    flag_mea=0; #determine when a measurement should take place
    E=np.double(0.0);
    W=np.double(0.0);

    #print('dtype O,w',O.dtype,w.dtype)
    #print('O',O)
    #print('w',w)

    iblk=0;
    
    
    # Preallocate arrays:
    E_blk=np.zeros((N_blk,1),dtype=np.double); # array to store the energy measured in every block
    W_blk=np.zeros((N_blk,1),dtype=np.double); # array to store the total weight in every block

    ## Equilibration phase
    for i_blk in range(0, N_eqblk):
        for j_step in range(0,N_blksteps):
            [Phi, w, O, E, W] = stpw.stepwlk(Phi, N_wlk, N_sites, w, O, E, W, H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld);
            if np.mod(j_step,itv_modsvd)==0:
                [Phi, O] = stb.stblz(Phi, N_wlk, O, N_up, N_par); # re-orthonormalize the walkers
            if np.mod(j_step,itv_pc)==0:
                [Phi, w, O]=pc.pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par); # population control


    ## Measurement phase    
    for i_blk in range(0,N_blk):
        for j_step in range(0,N_blksteps):
            if np.mod(j_step,itv_Em)==0:
                flag_mea=1;
            else:
                flag_mea=0;
                
            # propagate the walkers:
            [Phi, w, O, E_blk[i_blk], W_blk[i_blk]] = stpw.stepwlk(Phi, N_wlk, N_sites, w, O, E_blk[i_blk], W_blk[i_blk], H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld);
            if np.mod(j_step,itv_modsvd)==0:
                [Phi, O] = stb.stblz(Phi, N_wlk, O, N_up, N_par); # re-orthonormalize the walkers

            if np.mod(j_step,itv_pc)==0:
                [Phi, w, O]=pc.pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par); # population control
        
            if np.mod(j_step, itv_Em)==0:
                # update the exponent of the pre-factor exp(-deltau*(H-E_T))
                fac_norm=(np.real(E_blk[i_blk]/W_blk[i_blk])-0.5*U*N_par)*deltau;
        
    
        E_blk[i_blk]=(1.0)*E_blk[i_blk]/W_blk[i_blk];
        #print("E_blk",E_blk[i_blk])
        print('E('+str(i_blk)+')='+str(np.real(E_blk[i_blk])));


    ## Results
    E=np.real(E_blk);
    E_ave=np.mean(E)
    E_err=np.std(E)/np.sqrt(N_blk)
    # The total computational time:
    finish_time = tm.time(); # stop the  timer
    time =finish_time-start_time

    ## Save data to a *.mat file
    #save (savedFileName, 'E', 'E_ave', 'E_err', 'time');
    #save (savedFileName, '-append', 'Lx', 'Ly','Lz', 'N_up', 'N_dn', 'kx', 'ky','kz', 'U', 'tx', 'ty','tz');
    #save (savedFileName, '-append', 'deltau', 'N_wlk', 'N_blksteps', 'N_eqblk', 'N_blk', 'itv_pc','itv_modsvd','itv_Em');
    #save (savedFileName, '-append', 'H_k', 'E_nonint_v', 'Phi_T');
    np.savez(savedFileName, E=E, E_ave=E_ave, E_err=E_err, time=time, Lx=Lx,Ly=Ly,Lz=Lz,N_up=N_up,N_dn=N_dn, kx=kx, ky=ky, kz=kz,U=U,tx=tx,ty=ty,tz=tz,deltau=deltau,N_wlk=N_wlk,N_blksteps=N_blksteps,N_eqblk=N_eqblk,N_blk=N_blk,itv_pc=itv_pc,itv_modsvd=itv_modsvd,itv_Em=itv_Em,H_k=H_k,E_nonint_v=E_nonint_v,Phi_T=Phi_T); 
    
    

    ## Explanation of saved quantities:
    # E: the array of energy of each block
    # time: The total computational time
    # E_nonint_v: the non-interacting energy levels of the system
    # Phi_T: the trial wave function
    # For other saved quantities, type "help CPMC_Lab"

    return E_ave, E_err, savedFileName

