# Script to initialize internal quantities

import numpy as np
import scipy.linalg as spalg
import time as tm
import H_K as obh
import validation as valid


def initialization(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix):
    ## Check the validity of user inputs
    valid.validation(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix);


    ## Initialize internal quantities
    N_sites=Lx*Ly*Lz;
    N_par=N_up+N_dn;
    # form the one-body kinetic Hamiltonian
    if (np.any([np.absolute(kx),np.absolute(ky),np.absolute(kz)]))>0.0:
        H_k=np.zeros((N_sites,N_sites),dtype=np.cdouble);
    else:
        H_k=np.zeros((N_sites,N_sites),dtype=np.double);
    #H_k=np.zeros((N_sites,N_sites));
    obh.H_K(H_k,Lx, Ly,Lz, kx, ky,kz, tx, ty,tz);
    #print("H_k",H_k)
    # the matrix of the operator exp(-deltau*K/2)
    Proj_k_half = spalg.expm(-0.5*deltau*H_k); 
    #print("Proj_k_half",Proj_k_half)
    ## Initialize the trial wave function and calculate the ensemble's initial energy 
    # Diagonalize the one-body kinetic Hamiltonian to get the non-interacting single-particle orbitals:
    [E_nonint_v,psi_nonint] = np.linalg.eigh(H_k);
    #print("psi_noint",psi_nonint.dtype)
    # assemble the non-interacting single-particle orbitals into a Slater determinant:
    Phi_T=np.array(np.hstack((psi_nonint[:,:N_up],psi_nonint[:,:N_dn])),dtype=H_k.dtype);
    #print("Phi_T type",Phi_T.dtype)
    #print('Phi_T',Phi_T.conj().T)
    # the kinetic energy of the trial wave function
    E_K=np.sum(E_nonint_v[:N_up])+np.sum(E_nonint_v[:N_dn]);
    # the potential energy of the trial wave function
    n_r_up=np.diag(Phi_T[:,:N_up].dot(Phi_T[:,:N_up].conj().T));
    n_r_dn=np.diag(Phi_T[:,N_up:N_par].dot(Phi_T[:,N_up:N_par].conj().T));
    E_V=U*(n_r_up.conj().T.dot(n_r_dn));
    # the total energy of the trial wave function = the initial trial energy
    E_T = E_K+E_V;
    #print("E_T",E_T);
    
    ## Assemble the initial population of walkers
    #Phi=np.zeros((N_sites,N_par,N_wlk), dtype=np.cdouble);
    Phi=np.zeros((N_sites,N_par,N_wlk), dtype=H_k.dtype);

    # initiate each walker to be the trial wave function
    for i in range(0,N_wlk):
        # Phi(:,:,i) is the ith walker. Each is a matrix of size N_sites by N_par
        # The left N_sites by N_up block is the spin up sector
        # The rest is the spin down sector
        # They are propagated independently and only share the auxiliary field
        Phi[:,:,i]=Phi_T; 

    # initiate the weight and overlap of each walker to 1
    w=np.ones((N_wlk,1),dtype=np.double);
    #O=np.ones((N_wlk,1),dtype=np.cdouble);
    O=np.ones((N_wlk,1),dtype=H_k.dtype);
    # the arrays that store the energy and weight at each block
    E_blk=np.zeros((N_blk,1),dtype=np.double);
    W_blk=np.zeros((N_blk,1),dtype=np.double);

    ## initialize auxiliary filed constants
    # exponent of the prefactor exp(-deltau*(-E_T)) in the ground state projector 
    # fac_norm also include -0.5*U*(N_up+N_dn), the exponent of the prefactor in the Hirsch transformation
    fac_norm=(np.real(E_T)-0.5*U*N_par)*deltau; 
    #gamma in Hirsch's transformation
    gamma=np.arccosh(np.exp(0.5*deltau*U)); 
    # aux_fld is the 2x2 matrix containing all the possible values of the quantity exp(-gamma*s(sigma)*x_i)
    aux_fld=np.zeros((2,2), dtype=np.double); 
    # The first index corresponds to spin up or down
    # The second index corresponds to the auxiliary field x_i=1 or x_i=-1
    for i in range(0,2):
        for j in range(0,2):
            aux_fld[i,j]=np.exp(gamma*np.power(-1,i+j));
    #print("aux_fld",aux_fld)
    ## filename to be saved
    savedFileName=str(Lx)+'x'+str(Ly)+'x'+str(Lz)+'_'+str(N_up)+'u'+str(N_dn)+'d_U'+str(U)+'_kx'+str(kx)+'_ky'+str(ky)+'_kz'+str(kz)+'_Nwlk_'+str(N_wlk)+suffix+'.py';



    ## randomize the random number generator seed based on the current time
    np.random.seed(seed=np.int64(tm.time()));


    #returned variables
    #N_sites,N_par,H_k,Proj_k_half,E_nonint_v,Phi_T,E_T, Phi,W,O,fac_norm,aux_fld,savedFileName 
    
    return N_sites,N_par,H_k,Proj_k_half,E_nonint_v,Phi_T,E_T, Phi,w,O,fac_norm,aux_fld,savedFileName 
