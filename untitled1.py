#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 20:07:18 2017

@author: conor
"""

import numpy as np
import matplotlib.pyplot as plt
import mypbm

#2.26918531421
J = 1.0 # interaction energy

m = n = 100 # matrix dimensions

T =  2.269 # temperature
kb = 1.0 #1.38e-23 # boltzmann constant
B = 1.0/(kb*T) # Beta

def draw_mat(m,n): # to create a matrix
    z = np.random.random_sample(m*n).reshape(m,n)
    for i in range(n):
        for j in range(n):
            if z[i,j] >= 0.5:
                z[i,j] = 1
            elif z[i,j] < 0.5:
                z[i,j] = -1
                
    return z

g = draw_mat(m,n) # initial matrix

#plt.matshow(g)

def hamilton(g): # to calculate delta H = energy change due to flip
    q = np.copy(g)
    for i in range(n):
        for j in range(n):
            if i != 0 and j != 0 and i != len(g)-1 and j != len(g)-1:
                q[i,j] = 2*g[i,j]*(g[i+1,j]+g[i-1,j]+g[i,j+1]+g[i,j-1])
            elif i == 0 and j == 0: # periodic boundary conditions
                q[i,j] = 2*g[i,j]*(g[i+1,j]+g[i,j+1]+g[i,len(g)-1]+g[len(g)-1,j])
            elif i == 0 and j == len(g)-1:
                q[i,j] = 2*g[i,j]*(g[i,j-1]+g[i+1,j]+g[len(g)-1,len(g)-1]+g[0,0])
            elif i == len(g)-1 and j == 0:
                q[i,j] = 2*g[i,j]*(g[i-1,j]+g[i,j+1]+g[len(g)-1,len(g)-1]+g[0,0])
            elif i == len(g)-1 and j == len(g)-1:
                q[i,j] = 2*g[i,j]*(g[i,j-1]+g[i-1,j]+g[len(g)-1,0]+g[0,len(g)-1])
            elif i == 0 and j != 0 and j != len(g)-1:
                q[i,j] = 2*g[i,j]*(g[i,j-1]+g[i,j+1]+g[i+1,j]+g[len(g)-1,j])
            elif i == len(g)-1 and j != 0 and j != len(g)-1:
                q[i,j] = 2*g[i,j]*(g[i,j-1]+g[i,j+1]+g[i-1,j]+g[0,j])
            elif i != 0 and i != len(g)-1 and j == 0:
                q[i,j] = 2*g[i,j]*(g[i-1,j]+g[i+1,j]+g[i,j+1]+g[i,len(g)-1])
            elif i != 0 and i != len(g)-1 and j == len(g)-1:
                q[i,j] = 2*g[i,j]*(g[i+1,j]+g[i-1,j]+g[i,j-1]+g[i,0]) 
    for i in range(n):
        for j in range(n):
            q[i,j] *= g[i,j]
    return q

q = hamilton(g) # matrix of delta H values

def flippin(q): # function to decide whether to flip spins
    s = np.copy(g)
    for i in range(n):
        for j in range(n):
            if q[i,j] <= 0 or np.exp(-q[i,j]/T) > np.random.random():
                s[i,j] *= -1
    return s

g = flippin(q) # updated matrix


for iter in range(1000): # to iterate the above procedures
    q = hamilton(g)
    g = flippin(q)
    


for i in range(n): # set to 1 and 0 for pbm file
    for j in range(n): # should be commented out when calculating 
        if g[i,j] == -1: # data below
            g[i,j] = 0



mypbm.myplot(g,"2_269matrix.pbm")

"""
To measure the mean of absolute value of magnetisation per 
site <|M|> at different temperatures
"""

t_arr = np.linspace(0.4,3.0,40)
m_arr = []
k=0
ener = 0
ener_arr = []
ener2 = 0
ener2_arr = []

for T in t_arr:                     # function to iterate over matrix and 
    for iter in range(100):         # calculate magnetisation and 
        q = hamilton(g)             # energy per matrix
        g = flippin(q)
    for i in range(n):
        for j in range(n):
            k += g[i,j]
            k = abs(k) 
            ener += (q[i,j])/(-2*g[i,j])
            ener2 += (q[i,j])/(-2*g[i,j])
    m_arr.append(k)
    ener_arr.append(ener)
    ener2_arr.append(ener2)
    k = 0
    ener = 0
    ener2 = 0

m_arr = [x/(n**2) for x in m_arr]


plt.scatter(t_arr,m_arr) # plot magnetisation per site vs temp
plt.axis([0.0,3.4,0.0,1.0])
plt.title('Magnetisation per Site vs Temperature')
plt.xlabel('Temperature [J/kb]')
plt.ylabel(r'|<M>| per site [$\mu$]')

#plt.matshow(g)

"""
Unsuccessful attempt to calculate specific heat
"""

#ener_arr = [x/4 for x in ener_arr]
#ener2_arr = [(x/4)**2 for x in ener2_arr]
#ener2_arr = [x/500 for x in ener2_arr]
#ener_arr = [x**2/(500**2) for x in ener_arr]
#delE = [x-y for x,y in zip(ener2_arr,ener_arr)]
#
#hc_arr = [x/((y**2)*(n**2)) for x,y in zip(delE,t_arr)]
            
#plt.scatter(t_arr,hc_arr)
            
""" 
To calculate the energy per site
"""

ener_arr = [x/(2*(n**2)) for x in ener_arr]

plt.figure()
plt.scatter(t_arr,ener_arr)     # plot energy per site vs temp
plt.title('Energy per Site')
plt.xlabel('Temperature [J/kb]')
plt.ylabel('<E> per site [J]')
     
# Ising_Model
