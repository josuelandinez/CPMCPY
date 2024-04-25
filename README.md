# CPMCPY: Constrained Path Monte Carlo in Python  

CPMCPY is an **Auxiliary Field Quantum Monte Carlo** (AFQMC) code for simulation of the Hubbard Model with Periodic and Twisted Boundary conditions based on CPMC-Lab:

homepage: https://cpmc-lab.wm.edu/ 
https://www.huy.dev/CPC.pdf

## Overview

It is a Constrained Path Aulxiliary Field Quantum Monte Carlo code in python for the Hubbard model. It could be used with any version of python>=2.7 and numpy.  

The source code have the same structure that CPMC-Lab ported to python, the output data is saved in a npz file (numpy binary compressed files) the total energy and other parameters of the simulation. 

The aim of this code is to provide a simple implementation of Auxiliary Field Qunautm Monte Carlo that allows the user to experiment with new algorithms and features to improve the accuracy and performance of the AFQMC methods. 

The code is currently under development new features and changes will be constinously added. 
