import numpy as np
import halfK as h2K
import V as Vp
import measure as ms
def stepwlk(phi, N_wlk, N_sites, w, O, E, W, H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld):
    # function [phi, w, O, E, W] = stepwlk(phi, N_wlk, N_sites, w, O, E, W, H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld)
    # Perform one step of the random walk
    # Inputs:
    #   phi: the whole ensemble of walkers
    #   N_wlk: the number of walkers
    #   N_sites: the total number of lattice sites
    #   w: the array of weights of all the walkers
    #   O: the array of overlaps of all the walkers
    #   E: the total energy of all walkers
    #   W: the total weight of all walkers
    #   H_k: the one-body kinetic Hamiltonian
    #   Proj_k_half: the matrix of the operator exp(-deltau*K/2)
    #   flag_mea: the flag (1 or 0) that specifies whether the energy should the measured in this step
    #   Phi_T: the matrix of the trial wave function
    #   N_up: the number of spin up electrons
    #   N_par: the total number of electrons
    #   U: the on-site repulsion strength in the Hubbard model
    #   fac_norm: the exponent of the pre-factor exp(-deltau*(H-E_T))
    #   aux_fld: the 2x2 matrix containing all the possible values of the quantity exp(gamma*s(sigma)*x_i) (used in V.m only)
    # Outputs:
    #   phi: the ensemble of walkers after propagation
    #   w: the new array of weights of all walkers
    #   O: the new array of overlaps of all walkers
    #   E: the new total energy of all walkers
    #   W: the new total weight of all walkers
    #   
    

    ## Propagate each walker:
    e=np.zeros((N_wlk,1),dtype=np.double); # Array containing the energy of each walker
    for i_wlk in range(0,N_wlk):
        #Phi=np.array(phi[:,:,i_wlk],dtype=complex);
        Phi=phi[:,:,i_wlk];
        if w[i_wlk]>0:
            # multiply by the pre-factor exp(-deltau*(E_T)) in the ground-state projector 
            # and by the prefactor exp(-0.5*U*(N_up+N_dn)) in the Hirsch transformation
            w[i_wlk]=w[i_wlk]*np.exp(fac_norm);
            # propagate by the kinetic term exp(-1/2*deltau*K)
            [Phi, w[i_wlk], O[i_wlk], invO_matrix_up, invO_matrix_dn]=h2K.halfK(Phi, w[i_wlk], O[i_wlk], Proj_k_half, Phi_T, N_up, N_par);
            if w[i_wlk]>0:
                # propagate each lattice site of a walker by the potential term:
                for j_site in range(0,N_sites):
                    if w[i_wlk]>0:
                        [Phi[j_site,:], O[i_wlk], w[i_wlk], invO_matrix_up, invO_matrix_dn]=Vp.V(Phi[j_site,:], Phi_T[j_site,:], N_up, N_par, O[i_wlk], w[i_wlk], invO_matrix_up, invO_matrix_dn, aux_fld);
            
            
            if w[i_wlk]>0:
                # propagate by the kinetic term exp(-1/2*deltau*K)
                [Phi[:,:], w[i_wlk], O[i_wlk], invO_matrix_up, invO_matrix_dn]=h2K.halfK(Phi[:,:], w[i_wlk], O[i_wlk], Proj_k_half, Phi_T, N_up, N_par);            
                if w[i_wlk]>0:
                    # measure the energy if needed:
                    if flag_mea==1:
                        e[i_wlk]=ms.measure(H_k, Phi[:,:], Phi_T,  invO_matrix_up, invO_matrix_dn, N_up, N_par, U);
                
        phi[:,:,i_wlk]=Phi;


    ## Compute the ensemble's total energy and weight if measurement took place
    if flag_mea==1:
        for i_wlk in range(0,N_wlk):
            if w[i_wlk]>0:
                E+=(e[i_wlk]*w[i_wlk]);
                W+=w[i_wlk];



    return phi, w, O, E, W
