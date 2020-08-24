import numpy as np

def halfK(phi, w, O, Proj_k_half, Phi_T, N_up, N_par):
# function [phi, w, O, invO_matrix_up, invO_matrix_dn] = halfK(phi, w, O, Proj_k_half, Phi_T, N_up, N_par)
# Propagate a walker by the kinetic energy propagator exp(-deltau*K/2)
# Inputs:
#   phi: the matrix of a single walker
#   w: the weight of that walker
#   O: the overlap of that walker
#   Proj_k_half: the matrix of the operator exp(-deltau*K/2)
#   Phi_T: the matrix of the trial wave function
#   N_up: the number of spin up electrons
#   N_par: the total number of electrons
# Outputs:
#   phi: the matrix of the propagated single walker
#   w: the weight of the propagated walker
#   O: the overlap of the propagated walker   
#   invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
#   invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix 

    ## propagate the walker by exp(-deltau*K/2)
    phi=Proj_k_half.dot(phi);

    ## update the inverse of the overlap
    invO_matrix_up=np.linalg.inv(Phi_T[:,:N_up].conj().T.dot(phi[:,:N_up]));
    #print("invO_matrix_up=",invO_matrix_up)
    invO_matrix_dn=np.linalg.inv(Phi_T[:,N_up:N_par].conj().T.dot(phi[:,N_up:N_par]));
    #print("invO_matrix_up=",invO_matrix_up)
    # calculate the new overlap
    O_new=np.real(1.0/(np.linalg.det(invO_matrix_up)*np.linalg.det(invO_matrix_dn)));
    #print("O_new=",O_new)
    O_ratio=np.real((1.0*O_new)/O);
    #print("Onew,Oratio",O_new,O_ratio,'\n')

    
    ## enforce the constrained path condition
    # If the new weight is negative (O_raio<0), kill the walker by setting its weight to zero
    # real(O_ratio) enforces the phase-free approximation in case of complex phases (because the condition O_ratio>0 only checks the real part of O_ratio)
    if O_ratio>0:
        O=O_new;
        w=w*np.real(O_ratio);
    else:
        w=0; 

    #print("types", phi.dtype,w.dtype,O.dtype,invO_matrix_up.dtype, invO_matrix_dn.dtype)    


    #print('halfk O',O.dtype,O)
    #print('halfk w',w.dtype,w)
    #print('halfk invO matrrix up',invO_matrix_up.dtype,invO_matrix_up)
    #print('halfk invO matrix dn',invO_matrix_up.dtype,invO_matrix_dn.dtype)
    
    return phi, w, O, invO_matrix_up, invO_matrix_dn
