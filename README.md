# CPMCPY 

#AFQMC code for the hubbard model in any dimension with Periodic and Twisted Boundary conditions based on CPMC-Lab:

homepage: [huynguyen.io/CPMC-Lab](https://www.huynguyen.io/CPMC-Lab).



## Overview

It is a very raw attempt it need to be structured using clases and objects in python, however it is very easy to modify all important objets are defined globally at the initialization then passed to the functions locally. It could be used with any version of python>=2.7 and numpy.  

The source code have the same structure that CPMC-Lab adapted to python, the npz are numpy binary compressed files to save the data the total energy and other parameters of the simulation. 