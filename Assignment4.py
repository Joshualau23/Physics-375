# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 15:47:50 2019

@author: joshu
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from math import pi,sqrt,log
import numpy as np
import scipy as sp
from scipy import constants as con
import matplotlib.pyplot as plt


def model(x,r):
    y = x[0]
    dy = x[1]
    xdot = [[],[],[]]
    xdot[0] = dy
    xdot[1] = - (4*pi*y**2 + ((2.0/r) - (1.0/y) * dy)*dy)
    dm = 4*pi*r**2*y 
    return [xdot[0],xdot[1]]

  
r = np.linspace(0.00000001,5,10000)
LE0 = odeint(model,[1,0],r)
rho_list = LE0[:,0]



fig = plt.figure()
plt.plot(r,rho_list,color='black')
#plt.ylim(-1,1)
plt.title("Density vs radius" )
plt.ylabel("Density (dimensionless)")
plt.xlabel('Radius (Dimensionless)')
plt.grid()
fig.savefig('2b1.png')
plt.show()


M0 = rho0[1]
fig = plt.figure()
plt.plot(r,M0,color='black')
#plt.ylim(-1,1)
plt.title("Mass vs radius" )
plt.ylabel("Mass (dimensionless)")
plt.xlabel('Radius (Dimensionless)')
plt.grid()
fig.savefig('2b2.png')
plt.show()


surfpress = []
radius = []
Msun = 1.9891e30
mu = 2.4
G = 6.67408e-11
mp = 1.6726219e-27
k = 1.38064852e-23
T = 10
AU = 1.496e11
for i in range(0,len(r)):
    a = ((M0[i]**2 * rho_list[i]) / Msun**2)*(1.0/G**3)*(k*T / (mu*mp))**4
    b = (r[i] * (G*Msun*mu*mp) / (k*T*M0[i])) / AU
    surfpress.append(a)
    radius.append(b)

fig = plt.figure()
plt.plot(radius,surfpress,color='black')
plt.xlim(0,70000)
plt.title("Surface Pressure vs radius" )
plt.ylabel("Surface Pressure (Pascals)")
plt.xlabel('Radius (AU)')
plt.grid()
fig.savefig('2c.png')
plt.show()


Mj = []



for i in range(1,len(r)):
    #m = 0.2*Msun*((rho_list*(Msun/M0[i])**2 *(k*T/ (G*mu*mp))**3 * (1/(3.0e-15)))**(-1/2.0))
    m = 0.2*1.989*10**30*((rho_list[i]*(1.989*10**30/M0[i])*(con.k*10/(con.G*1.989*10**15*2.4*con.m_p))**3)/(3.0*10**15))**(-1.0/2.0)
    Mj.append(m)

fig = plt.figure()
plt.plot(radius[1:],Mj,color='black')
plt.xlim(0,50000)
plt.title("Jean's Mass vs radius" )
plt.ylabel("Jeans Mass (kg)")
plt.xlabel('Radius (AU)')
plt.grid()
fig.savefig('2d.png')
plt.show()





