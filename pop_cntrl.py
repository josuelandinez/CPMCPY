import numpy as np

def pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par):
    # function [Phi, w, O]=pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par)
    # Perform population control with a simple "combing" method
    # Inputs:
    #   Phi: the whole ensemble of walkers
    #   w: array containing the weights of all the walkers
    #   O: array containing the overlaps of all the walkers
    #   N_wlk: the number of walkers
    #   N_sites: the total number of lattice sites
    #   N_par: the total number of electrons
    # Outputs:
    #   Phi: the new ensemble of walkers after population control
    #   w: the new array of weights
    #   O: the new array of overlaps
    #
   

    ## Preparation
    # Create empty matrices that will hold the outputs
    new_Phi=np.zeros((N_sites, N_par, N_wlk),dtype=np.cdouble); #in the end the number of walkers will still be N_wlk
    new_O=np.zeros((N_wlk,1),dtype=np.cdouble);
    # scaling factor to bring the current total weight back to the original level (=N_wlk)
    d=(1.0*N_wlk)/np.sum(w);
    # start the "comb" at a random position to avoid bias against the first walker
    sum_w=-1.0*np.random.rand();
    n_wlk=0;

    ## Apply the comb
    for i_wlk in range(0,N_wlk):
        sum_w+=(w[i_wlk]*d);
        #print("sum_w",sum_w)
        n=np.int32(np.squeeze(np.ceil(sum_w)));
        #print("n",n)
        for j in range(n_wlk,n):
            new_Phi[:,:,j]=Phi[:,:,i_wlk];
            new_O[j]=O[i_wlk];
            #print("j",j)
        n_wlk=n;


    ## Return the new population, weights and overlaps:
    Phi=new_Phi;
    O=new_O;
    # All new walkers have weights to 1 and the total weight = N_wlk
    w=np.ones((N_wlk,1),dtype=np.double);
    
    return Phi, w, O
    
