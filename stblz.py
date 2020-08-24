import numpy as np

def stblz(Phi, N_wlk, O, N_up, N_par):
    # function [Phi, O] = stblz(Phi, N_wlk, O, N_up, N_par)
    # Perform the modified Gram-Schmidt orthogonalization to stabilize the walkers
    # Inputs:
    #   Phi: the whole ensemble of walkers
    #   N_wlk: the number of walkkers
    #   O: the array of overlaps of all walkers
    #   N_up: the number of spin up electrons
    #   N_par: the total number of electrons
    # Outputs:
    #   Phi: the stabilized ensemble of walkers
    #   O: the updated array of overlaps

    ## Perform the QR decomposition on each walker
    # Keep only the Q matrices and discard the R matrices
    for i_wlk in range(0,N_wlk):
        # for the spin up sector:
        [Phi[:,:N_up,i_wlk],R_up]=np.linalg.qr(Phi[:,:N_up,i_wlk],mode='reduced');
    
        # for the spin down sector:
        [Phi[:,N_up:N_par,i_wlk],R_dn]=np.linalg.qr(Phi[:,N_up:N_par,i_wlk],mode='reduced');
        #print("R_dn type", R_dn.dtype,R_dn)
        # Update the weight of each walker
        O[i_wlk]=O[i_wlk]/np.linalg.det(R_up.real)/np.linalg.det(R_dn.real);
    
    return Phi, O
