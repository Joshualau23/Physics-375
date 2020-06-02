# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 21:42:36 2019

@author: joshu
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib
from math import pi,sqrt,log10
import numpy as np
import scipy as sp
from scipy import constants as con
import matplotlib.pyplot as plt


Rsun = 695510
Msun = 1.9891e30
hbar = con.hbar
G = con.G
Me = con.m_e
Mp = con.m_p
Mn = con.m_n

def R_wd(M):
    return (hbar**2 / (1000*G*Me*Mp**2))*(M*Msun/Mp)**(-1/3.0)
    #return (0.01*Rsun*(M/0.7)**(-1.0/3.0))

def R_ns(M):
    return (hbar**2 / (1000*G*Mn*Mp**2))*(M*Msun/Mp)**(-1/3.0)
    #return (11*(M/1.4)**(-1.0/3.0))

def R_bh(M):
    return 3*(M) / 1000

def R_density(M):
    rho = 1.5E5
    return (3*M*Msun /(4*pi*rho)**(1/3.0)) / 1000
    #return ((3*M/(4*np.pi*(0.599)))**(1.0/3.0)*0.02*Rsun)

M = np.linspace(0.01,5,1000)
logM = []
Rwhite = []
Rneutron = []
Rblack = []
Rdensity = []
for i in range(0,len(M)):
    a = log10(R_wd(M[i]))
    b = log10(R_ns(M[i]))
    c = log10(M[i])
    d = log10(R_bh(M[i]))
    f = log10(R_density(M[i]))
    Rwhite.append(a)
    Rneutron.append(b)
    Rblack.append(d)
    Rdensity.append(f)
    logM.append(c)


fig = plt.figure()
plt.plot(logM,Rwhite,color = "gold")
plt.plot(logM,Rneutron)
plt.plot(logM,Rblack)
plt.title("Radius vs Mass" )
plt.ylabel("log(R/km)")
plt.xlabel('log(M/Msun)')
plt.grid()
plt.axvspan(0, 0.5, alpha=0.5, color='grey')
plt.show()

with open('test.txt', 'w') as f:
    for item in M:
        f.write("%s\n" % item)

####2b


fig = plt.figure()
plt.plot(logM,Rwhite)
plt.plot(logM,Rneutron)
plt.plot(logM,Rblack)
plt.plot(logM,Rdensity, color = "gold")
plt.title("Radius vs Mass" )
plt.ylabel("log(R/km)")
plt.xlabel('log(M/Msun)')
plt.grid()
#fig.savefig('2b.png')
#plt.show()


###2c

M_wd = np.linspace(0.01, 1.4, 1000)
M_ns = np.linspace(1.4, 3, 1000)
M_bh = np.linspace(3, 5, 1000)

Rwhite2 = []
Rneutron2 = []
Rblack2 = []
Rdensity2 = []

for i in range(0,len(M_wd)):
    a = log10(R_wd(M_wd[i]))
    b = log10(R_ns(M_ns[i]))
    d = log10(R_bh(M_bh[i]))
    Rwhite2.append(a)
    Rneutron2.append(b)
    Rblack2.append(d)

M_wd = map(np.log10, M_wd)
M_ns = map(np.log10, M_ns)
M_bh = map(np.log10, M_bh)
    



fig = plt.figure()
plt.plot(M_wd,Rwhite2)
plt.plot(M_ns,Rneutron2)
plt.plot(M_bh,Rblack2)
plt.title("Radius vs Mass" )
plt.ylabel("log(R/km)")
plt.xlabel('log(M/Msun)')
plt.grid()
#fig.savefig('2c.png')
#plt.show()
    
####2d


def w_wd(M):
    return 9.261e-6*M**(-1.0/3.0)

def w_ns(M):
    return 4.8214e-12*M**(-1.0/3.0)*Rsun**2

w_bh = 2.704e-11*696342**2

M_wd = np.linspace(0.01, 1.4, 1000)
M_ns = np.linspace(1.4, 3, 1000)
M_bh = 3


wwhite = []
wneutron = []
for i in range(0,len(M_wd)):
    a = w_wd(M_wd[i])
    b = w_ns(M_ns[i])
    wwhite.append(a)
    wneutron.append(b)
    

fig = plt.figure()
plt.loglog(M_wd,wwhite)
plt.loglog(M_ns,wneutron)
plt.loglog(M_bh,w_bh,marker='o',color = "gold")
plt.title("Revolution vs Mass log-log Plot" )
plt.ylabel("Revolution (revs/min)")
plt.xlabel('M/Msun')
plt.grid()
#fig.savefig('2d.png')
#plt.show()


h = con.h
G = con.G
c = con.c
kb = con.k
hbar = con.hbar

def Mevap(t):
    return (t * (c**4 * h) / (2560.0 * pi**2 * 4.0 * G**2))**(1/3.0)

tevap = 13.7e9 * 3.154e7
Mevaporation = Mevap(tevap)
#print Mevaporation

#print ((c**4 * h) / (2560.0 * pi**2 * 4.0 * G**2))**(1/3.0)


def temp(M):
    return (hbar*c**2)/ (8*pi*kb*G*M)

temperature = temp(Mevap(tevap))
#print temperature





