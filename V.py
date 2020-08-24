import numpy as np

def V(phi, phi_T, N_up, N_par, O, w, invO_matrix_up, invO_matrix_dn, aux_fld):
    # function [phi, O, w, invO_matrix_up, invO_matrix_dn] = V(phi, phi_T, N_up, N_par, O, w, invO_matrix_up, invO_matrix_dn, aux_fld)
    # Sample the auxiliary field over a single lattice site for a single walker and propagate that walker by the potential energy propagator exp(-deltau*V)
    # Inputs:
    #   phi: a single row in the matrix of a walker (corresponding to the amplitude of all electrons over a single lattice site in that walker)
    #   phi_T: the matrix of the trial wave function
    #   N_up: the number of spin up electrons
    #   N_par: the total number of electrons
    #   O: the overlap of the aforementioned walker
    #   w: the weight of the aforementioned walker
    #   invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
    #   invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix 
    #   aux_fld: the 2x2 matrix containing all the possible values of the quantity exp(gamma*s(sigma)*x_i) (used in V.m only)
    # Outputs:
    #   phi: the propagated row of the aforementioned walker
    #   O: the overlap after propagation of the aforementioned walker
    #   w: the weight after propagation of the aforementioned walker
    #   invO_matrix_up: the updated inverse of the spin up sector of the walker's overlap matrix 
    #   invO_matrix_dn: the updated inverse of the spin down sector of the walker's overlap matrix 
    

    ## Pre-allocate matrices:
    Gii=np.zeros((2,1),dtype=np.cdouble);
    RR=np.ones((2,2),dtype=np.cdouble);
    #matone=RR;
    matone=np.ones((2,2),dtype=np.double);
    
    ## Calculate the Green's function
    temp1_up=np.array(phi[:N_up].dot(invO_matrix_up),dtype=np.cdouble);
    #print("temp1_up",temp1_up)
    temp1_dn=np.array(phi[N_up:N_par].dot(invO_matrix_dn),dtype=np.cdouble);
    #print("temp1_dn",temp1_dn)
    temp2_up=np.array(invO_matrix_up.dot(phi_T[:N_up].conj().T),dtype=np.cdouble);
    #print("temp2_up",temp2_up)
    temp2_dn=np.array(invO_matrix_dn.dot(phi_T[N_up:N_par].conj().T),dtype=np.cdouble);
    #print("temp2_dn",temp2_dn)
    Gii[0]=temp1_up.dot(phi_T[:N_up].conj().T);
    Gii[1]=temp1_dn.dot(phi_T[N_up:N_par].conj().T);
    RR=np.multiply((aux_fld-matone),np.hstack((Gii,Gii)))+matone;
    #print("RR", RR.shape)
    #print("RR", RR)
    ## Perform the importance sampling and propagate the walker
    # compute overlaps
    O_ratio_temp = np.multiply(RR[0,:],RR[1,:]);
        
    #print "O_ratio_temp",O_ratio_temp
    # 
    O_ratio_temp_real=np.double(np.maximum(np.real(O_ratio_temp),np.zeros((1,2),dtype=np.double)));
    #print "O_ratio_temp real", np.maximum(np.real(O_ratio_temp),np.zeros((1,2),dtype=np.double));
    #print 'shape ortempreal',O_ratio_temp_real.shape
    #print 'O rario temp shape', O_ratio_temp.shape
    # the normalization for the importance-sampled pdf
    sum_O_ratio_temp_real=np.double(O_ratio_temp_real[0,0]+O_ratio_temp_real[0,1]);
    # if both auxiliary fields lead to negative overlap then kill the walker
    if sum_O_ratio_temp_real<=0.0:
        w=0;
    
    if w>0.0:
        # Otherwise update the weight
        w=w*0.5*sum_O_ratio_temp_real;    
    
        if (1.0*O_ratio_temp_real[0,0])/sum_O_ratio_temp_real>=np.random.rand():
            x_spin=np.int32(0);
        else:
            x_spin=np.int32(1);
        
        # propagates the walker with the chosen auxiliary field
        phi[:N_up]=phi[:N_up]*aux_fld[0,x_spin];
        phi[N_up:N_par]=phi[N_up:N_par]*aux_fld[1,x_spin];    
        # Update the overlap using Sherman-Morrison
        O=O*np.real(O_ratio_temp[x_spin]);
        invO_matrix_up=invO_matrix_up+((1.0-aux_fld[0,x_spin])/RR[0,x_spin])*np.dot(temp2_up,temp1_up);
        invO_matrix_dn=invO_matrix_dn+((1.0-aux_fld[1,x_spin])/RR[1,x_spin])*np.dot(temp2_dn,temp1_dn);
        
    return phi, O, w, invO_matrix_up, invO_matrix_dn
