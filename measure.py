import numpy as np

def measure(H_k, phi, Phi_T, invO_matrix_up, invO_matrix_dn, N_up, N_par, U):
# function e = measure(H_k, phi, Phi_T, invO_matrix_up, invO_matrix_dn, N_up, N_par, U)
# Calculate the mixed estimator for the ground state energy of a walker
# Inputs:
#   H_k: the one-body kinetic Hamiltonian
#   phi: the matrix of a single walker
#   Phi_T: the matrix of the trial wave function
#   invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
#   invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix 
#   N_up: the number of spin up electrons
#   N_par: the total number of electrons of both spins
#   U: the on-site repulsion strength in the Hubbard model
# Outputs:
#   e: the mixed estimator for the ground state energy of the input walker

    ##  calculate the single-particle Green's function matrix for each spin:
    temp_up=phi[:,:N_up].dot(invO_matrix_up);
    temp_dn=phi[:,N_up:N_par].dot(invO_matrix_dn);
    G_up=temp_up.dot(Phi_T[:,:N_up].conj().T);
    G_dn=temp_dn.dot(Phi_T[:,N_up:N_par].conj().T);

    ## calculate the potential energy:
    n_int=np.diag(G_up).T.dot(np.diag(G_dn));
    potentialEnergy=n_int*U;

    ## calculate the kinetic energy:
    kineticEnergy=np.sum(np.multiply(H_k.T,G_up+G_dn)); 
    #print("kineticEnergy")
    #print(kineticEnergy)

    #print("potentialEnergy")
    #print(potentialEnergy)
    
    ## calculate the total energy:
    e=potentialEnergy+kineticEnergy;
    #print("e")
    #print(e)

    
    return np.real(e)
