#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 12:49:05 2024

@author: randallclark
"""
"""
The goal of this script is to solve the Burgers Equation using a Fourier Galerkin approach to the Spectral Methods of PDE solving
∂u/∂t + u∂u/∂x - v∂^2u/∂x^2 = 0

For a periodic bundary condition from 0 to 2pi, we represent u as a fourier expansion:
    u_N = SUM_k=(from -N to +N) {u_k(t)*exp(ikx)}

By multiplying by e^(-ixk) and integrating the Burgers equation, we can pickout the kth u value (after plugging in our expansion) and by the property
of orthogonallity, we can simplify this all down to a simple equation of ∂u_k/∂t = du_k/dt

du_k/dt = -k^2*v*u_k - i/(2pi) * convolution sum_(j+l=k) {u_j*u_l*l}      Equation (1)

This can be solved with ODEINT and we get our initial condition by again making use of the orthogonallity property
int(u_N(x,t=0)*exp(-ikx))dx = int( sum_j(u_k(t))*exp(+ijx) *exp(-ikx))dx = 2pi*u_k(t=0)
u_k(t=0) = 1/2pi * int(u_N(x,t=0)*exp(-ikx))dx 
if we choose u_n(x,t=0) to have an initial condition equalling sin(x/2) (which satisfies our periodic boundary condition), we get:
u_k(t=0) = 2/(pi*(1-4k^2))

Combine these two steps, add up the terms of u_k, and we can recreate u_N, our spectral approximation of u, at any time!
"""

import numpy as np
from scipy.integrate import odeint


class PDE_Burgers():
    '''
    Generate Burgers data on convection/diffusion flow
    '''
    def __init__(self, N=10,v=1):
        """Initialize the PDE_Burgers object
        Args:
            N (int)         : Cuttoff number for the spectral expansion order (k will range from -N to +N)
            v (float)       : The diffusion term
        """
        self.N = N
        self.v = v
    
    
    """
    Initial Condition for u_k
    So far, only a sin(x/2) initial condition has been implemented
    """
    def GenerateInitialCondition(self):
        #Here we choose an initial condition function and evaluate the value of u_k(t=0) for each 2N+1 k's
        #Modification: Split the real and imaginary into two different dimensions (4*N+2 total dimensions now)
        #ex: Re(u1),Im(u1),Re(u2),Im(u2),...,Im(u2N+1)
        
        Ini = np.zeros((2*(2*self.N+1)))
        for i in range(2*self.N+1):
            Ini[i*2] = self.sin_half_x_initialcondition(-self.N+i) #Ini values start at -N and go to +N
        
        return Ini

    def sin_half_x_initialcondition(self,k):
        #The initial condition for the fourier coefficients evaluates to a very simple function when the
        #initial u(t=0,x) = sin(x/2)
        return (2)/(np.pi*(1-4*k**2))
    
    """
    Integration tools for u_k
    """
    def f_spectral_burger(self,u,t):
        dukdt = []
        for k in range(-self.N,self.N+1):
            dudt = self.du_kdt(k,u)
            dukdt.append(np.real(dudt))
            dukdt.append(np.imag(dudt))
        return dukdt
    
    def du_kdt(self,k,u):
        #First step, include the first term in Equation 1 (Note the splitting of real and imaginary terms)
        Answer = -k**2*self.v*(u[2*(k+self.N)]+1j*u[2*(k+self.N)+1])
        
        #Now Sum over the terms in the convolutional sum
        Sum = 0
        j = -self.N
        l = self.N
        
        #Discount any terms that would exceed N in the convolutional sum
        if k > 0:
            j = j+k
        if k < 0:
            l = l+k
        
        #Perform the summation
        for i in range(2*self.N+1-abs(k)):
            Sum += (u[2*(j+self.N)] + 1j*u[2*(j+self.N)+1]) * ((u[2*(l+self.N)] + 1j*u[2*(l+self.N)+1])) * l
            j+=1
            l+=-1
        return Answer - 1j/(2*np.pi)*Sum

    """
    #Now we need a function that can add up all of the values! to get our original u(x,t) back!  
    """
    def rebuild_uN(self,u,x):
        #We take each (Re(u)+i Im(u))*(cos(kx)+)
        uN = 0 + 0.0j
        for i in range(2*self.N+1):
            k = self.N-i
            uN += (u[2*i]+1j*u[2*i+1])*(np.cos(x*k)-1j*np.sin(x*k))
        return np.real(uN)
    
    def sample_uN(self,u,L):
        y = np.zeros(L)
        x = np.linspace(0,2*np.pi,L)
        for i in range(L):
            y[i] = self.rebuild_uN(u,x[i])
        return y
        
    """
    Call this function to get data!
    """
    def generate_uN(self,InitialCondition,TimeRange,L = 100):
        #Integrate the 2*(2N+1) dimensions (real and imaginary) of the spectral expansion of uN, then add them all back up together to get uN
        Solver = odeint(self.f_spectral_burger,InitialCondition,TimeRange)
        uN = np.zeros((Solver.shape[0],L))
        for i in range(Solver.shape[0]):
            uN[i] = self.sample_uN(Solver[i],L)
        return uN
    
    def generate_Solver(self,InitialCondition,TimeRange,L = 100):
        #Integrate the 2*(2N+1) dimensions (real and imaginary) of the spectral expansion of uN, then add them all back up together to get uN
        Solver = odeint(self.f_spectral_burger,InitialCondition,TimeRange)
        return Solver
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    