# Script to check the validity of all user inputs
import numpy as np


def validation(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix):


    if np.mod(Lx,1)!=0 or Lx<=0 :
        print('Error: Lx must be a positive integer!')
        #break
        raise SystemExit
    elif np.mod(Ly,1)!=0 or Ly<=0 :
        print('Error: Ly must be a positive integer!')
        raise SystemExit
    elif np.mod(Lz,1)!=0 or Lz<=0 :
        print('Error: Lz must be a positive integer!')
        raise SystemExit
    elif np.mod(N_up,1)!=0 or N_up<0 :
        print('Error: N_up must be a non-negative integer!')
        raise SystemExit
    elif np.mod(N_dn,1)!=0 or N_dn<0 :
        print('Error: N_dn must be a non-negative integer!')
        raise SystemExit
    elif N_up > ( 2*Lx*Ly*Lz) :
        print('Error : there are too many spin-up electrons on this lattice')
        raise SystemExit
    elif N_dn > ( 2*Lx*Ly*Lz) :
        print('Error: there are too many spin-down electrons on this lattice')
        raise SystemExit
    elif kx<=-1.0 or kx>1.0 :
        print('Error: kx must be in the interval (-1,1]!')
        raise SystemExit
    elif ky<=-1.0 or ky>1.0 :
        print('Error: ky must be in the interval (-1,1]!')
        raise SystemExit
    elif kz<=-1.0 or kz>1.0 :
        print('Error: kz must be in the interval (-1,1]!')
        raise SystemExit
    elif U<0.0 :
        print('Error : U must be non-negative in the repulsive Hubbard model!')
        raise SystemExit
    elif tx<0.0 :
        print('Error: tx must be non-negative!')
        raise SystemExit
    elif ty<0.0 :
        print('Error: ty must be non-negative!')
        raise SystemExit
    elif tz<0.0 :
        print('Error: tz must be non-negative')
        raise SystemExit
    elif deltau<0.0 :
        print('Error: deltau must be positive!')
        raise SystemExit
    elif np.mod(N_wlk,1)!=0 or N_wlk<=0 :
        print('Error: N_wlk must be a positive integer!')
        raise SystemExit
    elif np.mod(N_blksteps,1)!=0 or N_blksteps<=0 :
        print('Error: N_blksteps must be a positive integer!')
        raise SystemExit
    elif np.mod(N_eqblk,1)!=0 or N_eqblk<0 :
        print('Error: N_eqblk must be a non-negative integer!')
        raise SystemExit
    elif np.mod(N_blk,1)!=0 or N_blk<=0 :
        print('Error: N_blk must be a positive integer!')
        raise SystemExit
    elif np.mod(itv_modsvd,1)!=0 or itv_modsvd<=0 :
        print('Error: itv_modsvd must be a positive integer!')
        raise SystemExit
    elif np.mod(itv_pc,1)!=0 or itv_pc<=0 :
        print('Error: itv_pc must be a positive integer!')
        raise SystemExit
    elif np.mod(itv_Em,1)!=0 or itv_Em<=0 or itv_Em>N_blksteps :
        print('Error: itv_Em must be a positive integer no greater than N_blksteps!')
        raise SystemExit
    elif not(isinstance(suffix,str)): # llok how to chek if a variable is a string in python
        print ('suffix must be a string!')
        raise SystemExit
    else:
        if deltau>1:
            print('Warning: imaginary time step deltau>1 is too large!!!')
            print('Ctrl+c to stop running!')
    
        if N_wlk*N_blksteps*(N_eqblk+N_blk)*Lx*Ly*(N_up+N_dn)>1e11 :
            print('Warning: N_wlk*N_blksteps*(N_eqblk+N_blk)*Lx*Ly*(N_up+N_dn) > 1e11')
            print('This computation might take more than a day!')
            print('Ctrl+c to stop running!')
    
        if itv_modsvd>N_blksteps:
            print('Warning: itv_modsvd is greater than N_blksteps. The walkers will not be periodically re-orthogonalized.')
            print('Ctrl+c to stop running!')
    
        if itv_pc>N_blksteps:
            print('itv_pc is greater than N_blksteps. There will be no population control.')
            print('Ctrl+c to stop running!')
    

