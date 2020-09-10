#!/usr/bin/env python
# coding: utf-8


# This will be an emulation of the Fortran code.
# Rather than using FOR loops, I hope to use IF statements to populate the matrix.
# This method should be faster.


import numpy as np
import concurrent.futures
#import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy import linalg as la
import time
import gc
gc.enable()

print("*************************************************************************")

# input data
Nr, Nphi, Nz = 2, 4, 2                     # number of rings and angles
Rext, Rint = 1.0, 0.8                         # exterior and interior radii
Nsites = Nr*Nphi                              # number of sites in 2D
Nsites_3D = Nr*Nphi*Nz		  	      # number of sites in 3D
Nstat = 2*Nsites                              # number of states in 2D with spin
Nstat_3D = 2*Nsites_3D			      # number of states in 3D with spin
ts = 1.0                                      # energy unit = hbar**2/(2*m*Rext**2)
gamma = -0.2                                  # InAs geff=-14.9, meff=0.023, gamma=geff*meff/2
alphaR, betaD = 0.7, 0.3                      # (SOI param)/Rext in units
#emag = 1.0                                   # Cyclotron energy in units ts

# single particle state initializers
pos = np.zeros([Nsites, 4])                   # position of sites in polar and cartesian coord.
pos3D = np.zeros([Nsites_3D, 6])                  # position of sites in 3D cyl. and Cart.
qn = np.zeros([Nstat, 2])                     # quantum numbers (2D)
qn3D = np.zeros([Nstat_3D, 3])		      # quantum numbers (3D)
Ea = np.zeros(Nstat)                          # single particle energies (2D)
Ea3D = np.zeros(Nstat_3D)                      # single particle energies (3D)
#Hmat = np.zeros([Nstat,Nstat], dtype=complex) # single particle Hamiltonian
Psi = np.zeros([Nstat,Nstat], dtype=complex)  # single particle eigenvector matrix
Pari = np.zeros(Nstat)                        # parity of single particle states

# defining sigma matrices
def SigMat(type, s1, s2, phi):                    # Pauli Matrices
    # Sigma x
    if (type == 'x'):
        if s1 == s2:
            SigMat = 0.0 + 0.0j                   # s1=s2=1 or s1=s2=-1
        else:
            SigMat = 1.0 + 0.0j                   # s1=-s2=1 or -s1=s2=1

    # Sigma y
    if (type == 'y'):
        if s1 == s2:
            SigMat = 0.0 + 0.0j                   # s1=s2=1 or s1=s2=-1
        else:
            SigMat = 0.0 + s2*1.0j                # s1=-s2=1 or -s1=s2=1

    # Sigma z
    if (type == 'z'):
        if s1 == s2:
            SigMat = s2*(1.0 + 0.0j)              # s1=s2=1 or s1=s2=-1
        else:
            SigMat = s2*(0.0 + 0.0j)              # s1=-s2=1 or -s1=s2=1

    # Sigma r = Sigma x Cos(phi) + Sigma y Sin(phi)
    if (type == 'r'):
        if s1 == s2:
            SigMat = 0.0 + 0.0j                   # s1=s2=1 or s1=s2=-1
        else:
            SigMat = np.exp(s2*phi*1.0j)          # s1=-s2=1 or -s1=s2=1

    # Sigma p = - Sigma x Sin(phi) + Sigma y Cos(phi)
    if (type == 'p'):
        if s1 == s2:
            SigMat = 0.0 + 0.0j                   # s1=s2=1 or s1=s2=-1
        else:
            SigMat = s2*1.0j*np.exp(s2*phi*1.0j)  # s1=-s2=1 or -s1=s2=1

    # Identity matrix
    if (type == '1'):
        if s1 == s2:
            SigMat = 1.0 + 0.0j                   # s1=s2=1 or s1=s2=-1
        else:
            SigMat = 0.0 + 0.0j                   # s1=-s2=1 or -s1=s2=1

    return SigMat

print(f"Creating Lattice of {Nstat_3D} Sites...")
# 3D Lattice
k = 0
r = Rext
d_r = 0
d_z = 1/Nz
if (Nr>1): d_r = (Rext-Rint)/(Nr-1)
d_phi = 2*np.pi/Nphi
for iz in np.arange(1, Nz+1, 1, dtype=int):
    z = (iz-1)*d_z
    for ir in np.arange(1, Nr+1, 1, dtype=int):
        if (Nr>1): r = Rext-d_r*(ir-1)
        for iphi in np.arange(1, Nphi+1, 1, dtype=int):
            phi = d_phi*(iphi-1)
            k+=1
            pos3D[k-1,0] = r
            pos3D[k-1,1] = phi
            pos3D[k-1,2] = z
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            pos3D[k-1,3] = x
            pos3D[k-1,4] = y
            pos3D[k-1,5] = z
            if (r==0):
                break

print("Defining Quantum Numbers...")
# Quantum Numbers 3D
a = 0
for k in np.arange(1, Nsites+1, 1, dtype=int):
    for spin in np.arange(1,-2,-2, dtype=int):
        a+=1
        qn3D[a-1,0] = k
        qn3D[a-1,1] = spin

print("Defining the Hamiltonian...")
def Hmat(emag):
    # Single Particle Hamiltonian
    H = np.zeros([Nstat_3D,Nstat_3D], dtype=complex)

    for a1 in np.arange(1, Nstat_3D+1, 1, dtype=int):
        #Hmat in 3D (needs another for loop on the outside for z slice)
        k1 = int(qn3D[a1-1,0])
        spin1 = int(qn3D[a1-1,1])
        r1 = pos3D[k1-1,0]
        phi1 = pos3D[k1-1,1]
        z1 = pos3D[k1-1,2]
        n1 = int(round((Rext-r1)/d_r+1))
        j1 = int(round(phi1/d_phi)+1)
        int_z1 = int(round(z1/d_z)+1)
        if (d_r == 0): n1=1
        tphi = ts/d_phi**2*(Rext/r1)**2
        tz = 1/d_z**2
        E = 0
        #H[a1-1,a1-1] = H[a1-1,a1-1] + 2*tphi + (emag*r1/4)**2 + gamma*emag*spin1/2
        # Hamiltonian in 3D with planar electric field

        H[a1-1,a1-1] = H[a1-1,a1-1] + 2*tphi + 2*tz + (emag*r1/4)**2 + gamma*emag*spin1/2 + E*r*(np.cos(phi1) + np.sin(phi1))

        for a2 in np.arange(a1, Nstat_3D+1, 1, dtype=int):
            k2 = int(qn3D[a2-1,0])
            spin2 = int(qn3D[a2-1,1])
            r2 = pos3D[k2-1,0]
            phi2 = pos3D[k2-1,1]
            z2 = pos3D[k2-1,2]
            n2 = int(round((Rext-r2)/d_r+1))
            j2 = int(round(phi2/d_phi)+1)
            int_z2 = int(round(z2/d_z)+1)
            if (d_r == 0): n2=1

            # rr = abs(r1-r2)
            # pp = abs(phi1-phi2)
            nn = abs(n1-n2)
            jj = abs(j1-j2)
            #print(r1,r2)
            if (spin1==spin2):
                if (abs(int_z1-int_z2)==0):
                    if (n1==n2):
                        if (jj==1 or jj==(Nphi-1)): H[a1-1,a2-1]=H[a1-1,a2-1]-tphi
                    if (Nr>1):
                        tr = ts*(Rext/d_r)**2
                        if (jj==0 and nn==0): H[a1-1,a2-1]=H[a1-1,a2-1]+2*tr
                        if (jj==0 and nn==1): H[a1-1,a2-1]=H[a1-1,a2-1]-tr
	                   # print( H[a1-1,a1-1] )
                        if (n1==n2 and (jj==1 or jj==(Nphi-1))):
                            semn=1
                            if ((jj==1 and j1>j2) or (j1==1 and j2==Nphi)): semn=-1
                            if (spin1==spin2): H[a1-1,a2-1]=H[a1-1,a2-1]-0.25*1.0j*emag/d_phi*semn
                if (int_z2-int_z1==1):
                    H[a1-1,a2-1]=H[a1-1,a2-1]-2*tz
                #print(pos3D[k1-1,2], pos3D[k2-1,2])
	            #print( H[a1-1,a2-1] )
        H[a2-1,a1-1]=np.conjugate(H[a1-1,a2-1])
	            #print(H[a1-1,a2-1])
    return H

print("*************************************************************************")

print(Hmat(0))

emag = [i for i in np.arange(0, 10.025, 0.025)]

''' This is the code for getting eigensystem in serial. Can take up to 10x as long.

print(f"Evaluating, in serial, the Hamiltonian at B-field strengths {emag[0]} to {round(emag[-1])}...")
start = time.time()

H_list = [Hmat(i) for i in np.arange(0,10,0.025)]

end = time.time()
print(f" Finished in {round(end - start, 2)} seconds.")

print("Calculating, in serial, Eigensystem...")
start = time.time()

Ea = [la.lapack.zheevd(H_list[i]) for i in range(400)]

end = time.time()
print(f" Finished in {round(end - start, 2)} seconds.")

print("*************************************************************************")
'''

print(f"Evaluating, in parallel, the Hamiltonian at B-field strengths {emag[0]} to {round(emag[-1])}...")
start = time.time()
with concurrent.futures.ProcessPoolExecutor() as pool:
    out = pool.map(Hmat, emag)
end = time.time()

print(f" Finished in {round(end - start, 2)} seconds.")

gc.collect()

'''
# The Eigensystem calculation Hangs on evaluation. I suspect from using low-level functions.
print("Calculating, in parallel, Eigensystem...")
start = time.time()
lowest_twelve = [value for value in out]
with concurrent.futures.ThreadPoolExecutor() as pool:
    eigsys = pool.map( la.lapack.zheevd, lowest_twelve )
end = time.time()
print(f" Finished in {round(end - start, 2)} seconds.")

'''
print("Calculating, in serial, Eigensystem...")
start = time.time()
Ea = []
for value in out:
    Ea.append(la.lapack.zheevd(value))
end = time.time()
print(f" Finished in {round(end - start, 2)} seconds.")
print(f"The number of eigenvalues is: {Ea[0][0].shape}")
print("*************************************************************************")

print(f"{Ea[0][0]}")

''' Since parallel Eigsys evaluation hangs, I will move on
'''











