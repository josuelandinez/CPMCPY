import numpy as np
import scipy as sp

def H_K(H,Lx,Ly,Lz,kx,ky,kz,tx,ty,tz): 
# function H=H_K(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz)
# Generate the one-body kinetic term of the Hubbard Hamiltonian with the given parameters
# Input:
#   Lx: The number of lattice sites in the x direction.
#   Ly: The number of lattice sites in the y direction.
#   Lz: The number of lattice sites in the z direction.
#   kx: The x component of the twist angle in TABC (twist-averaging boundary conditions)
#   ky: The y component of the twist angle in TABC
#   kz: The z component of the twist angle in TABC
#   tx: The hopping amplitude between nearest-neighbor sites in the x direction
#   ty: The hopping amplitude between nearest neighbor sites in the y direction
#   tz: The hopping amplitude between nearest neighbor sites in the y direction
# Output
#   H: The one-body kinetic Hamiltonian in the form of a square matrix of size (Lx*Ly*Lz) 
 

    r=0; # Lattice Site counter
    N_sites=Lx*Ly*Lz;
    kx=np.pi*kx*1.j;
    ky=np.pi*ky*1.j;
    kz=np.pi*kz*1.j;


    for mz in range(0,Lz):
        for iy in range(0,Ly):
            for jx in range(0,Lx):
                if Lx!=1:
                    if jx==0:
                        H[r,r-1+Lx]+=-tx*np.exp(kx);
                        H[r,r+1]+=-tx;
                    elif jx==(Lx-1):
                        H[r,r-1]+=-tx;
                        H[r,r+1-Lx]+=-tx*np.exp(-kx);
                    else:
                        H[r,r-1]=-tx;
                        H[r,r+1]=-tx;
                if Ly!=1:
                    if iy==0:
                        H[r,r+(Ly-1)*Lx]+=-ty*np.exp(ky);
                        H[r,r+Lx]+=-ty;
                    elif iy==(Ly-1):
                        H[r,r-Lx]+=-ty;
                        H[r,r-(Ly-1)*Lx]+=-ty*np.exp(-ky);
                    else:
                        H[r,r-Lx]=-ty;
                        H[r,r+Lx]=-ty;
                if Lz!=1:
                    if mz==0:
                        H[r,r+(Lz-1)*Lx*Ly] +=-tz*np.exp(kz);
                        H[r,r+Lx*Ly]+=-tz;
                    elif mz==(Lz-1):
                        H[r,r-Lx*Ly]+=-tz;
                        H[r,r-(Lz-1)*Lx*Ly]+=-tz*np.exp(-kz);
                    else:
                        H[r,r-Lx*Ly]=-tz;
                        H[r,r+Lx*Ly]=-tz;
                        
                r=r+1;      

