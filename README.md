# PDE-Spectral Method-Fourier Galerkin Approach-Burgers' Equation
This project was performed to practice a common approach to using the Spectral Method with the fourier Galerkin approach to solving PDE's. Specifically, we'll be tackling the periodic Burgers' Equation:
https://en.wikipedia.org/wiki/Burgers%27_equation

The approach is based off of Chapter 10, section 10.3.1 Fourier Galerkin in the book Computational methods in plasma physics.
Jardin, Stephen. Computational methods in plasma physics. CRC press, 2010.

------------------------------------------------------------------------------------------------------------------------------------------------
The goal of this Project is to solve the Burgers Equation using a Fourier Galerkin approach to the Spectral Methods of PDE solving
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

https://github.com/user-attachments/assets/8841bf6c-bc72-42fc-80e1-58845b050bbf
