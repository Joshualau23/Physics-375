# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 15:10:31 2019

@author: joshu
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from math import pi,sqrt,log


def model(x,t,k):
    y = x[0]
    dy = x[1]
    xdot = [[],[]]
    xdot[0] = dy
    xdot[1] = (-2.0/t)*dy - y**k
    return xdot




time = np.linspace(.000001,10,1000)


k=0
LE0 = odeint(model,[1,0],time,args=(k,))

k=1
LE1 = odeint(model,[1,0],time,args=(k,))
k=2
LE2 = odeint(model,[1,0],time,args=(k,))
k=3
LE3 = odeint(model,[1,0],time,args=(k,))
k=4
LE4 = odeint(model,[1,0],time,args=(k,))
k=5
LE5 = odeint(model,[1,0],time,args=(k,))

fig = plt.figure()
plt.plot(time,LE0[:,0],color='black')
plt.plot(time,LE1[:,0],color='red')
plt.plot(time,LE2[:,0],color='yellow')
plt.plot(time,LE3[:,0],color='green')
plt.plot(time,LE4[:,0],color='blue')
plt.plot(time,LE5[:,0],color='purple')
plt.ylim(-1,1)
plt.title("$\Theta$ versus x" )
plt.xticks(np.arange(0, 10, 1))
plt.ylabel("$\Theta$")
plt.xlabel('$x$')
plt.grid()
#fig.savefig('3e.png')
plt.show()


plt.plot(time,LE3[:,0],color='Gold')
plt.plot(time,LE3[:,1],color='pink')
plt.ylim(-1,1)
plt.title("$\Theta$ versus x" )
plt.xticks(np.arange(0, 10, 1))
plt.ylabel("$\Theta$")
plt.xlabel('$x$')
plt.grid()
plt.show()


###g


time = np.array(time)
theta = np.array(LE3[:,0])
dtheta = np.array(LE3[:,1])

x0 = 0
dthetadx = 0
for i in range(0,len(time)):
    a = round(theta[i],8)
    if a == -1.98e-06:
        x0 = time[i]
        dthetadx = dtheta[i]
    else:
        continue
        

print x0
print dthetadx




Msun = 0.5*1.9891e30
Rsun = 696000000*0.6
rhoc = Msun / (-4*pi*(Rsun**3 / x0) * dthetadx)
print rhoc / 10**5

G = 6.67408e-11
n=3.0

K = (Rsun**2 / x0**2)*(4*pi*G) / ((1.0/n + 1)*n*rhoc**(1/n -1))
print K / 10**9


N = ((4*pi)**(1/3.0) / (n+1)) * (-x0**2 * dthetadx)**( (1-n) / n) * x0**((n-3) / n)
K1 =  Msun**(2/3.0)*N*G
print K1 / 10**9


###h

alpha = Rsun / x0
x = np.array(time)
rRstar = []
rho = []
for i in range(0,len(x)):
    a = x[i]*alpha / Rsun
    b = rhoc*theta[i]**3
    rRstar.append(a)
    rho.append(b)

fig = plt.figure()
plt.plot(rRstar,rho,color='royalblue')
plt.title("Density versus radius" )
#plt.xticks(np.arange(0, 10, 1))
plt.ylabel("Density")
plt.xlabel('r/Rsun')
plt.xlim(0,1)
plt.grid()
#fig.savefig('3h1.png')
plt.show()



P = []

for i in range(len(rRstar)):
    a = K*rho[i]**(1 + 1.0/n)
    P.append(a)

fig = plt.figure()
plt.plot(rRstar,P,color='lawngreen')
plt.title("Pressure versus radius" )
#plt.xticks(np.arange(0, 10, 1))
plt.ylabel("Pressure")
plt.xlabel('r/Rsun')
plt.xlim(0,1)
plt.grid()
#fig.savefig('3h2.png')
plt.show()

mu = 1/ (2*.55 + 0.75*0.4 + 0.5*0.05)
k = 1.38064852e-23
mp = 1.6726219e-27

T = [] 
for i in range(0,len(rho)):
    a = (mu*mp*P[i]) / (rho[i]*k) 
    T.append(a)

fig = plt.figure()
plt.plot(rRstar,T ,color='salmon')
plt.title("Temperature versus radius" )
#plt.xticks(np.arange(0, 10, 1))
plt.ylabel("Temperature")
plt.xlabel('r/Rsun')
plt.xlim(0,1)
plt.grid()
#fig.savefig('3h3.png')
plt.show()

print T[0] 

###i


X = 0.55
def epp(p,t):
    return 1.07e-7*X**2 * (p / 10**5)*(t/10**6)**4
energypp = []

for i in range(0,len(T)):
    a = epp(rho[i],T[i])
    energypp.append(a)

fig = plt.figure()
plt.plot(rRstar,energypp,color='salmon')
plt.title("Energy generation rate (p-p chain) versus radius" )
#plt.xticks(np.arange(0, 10, 1))
plt.ylabel("Energy generation rate")
plt.xlabel('r/Rsun')
plt.xlim(0,1)
plt.grid()
#fig.savefig('3i1.png')
plt.show()


Xcno = 0.03*X

def ecno(p,t):
    return 8.24e-26*X*Xcno * (p / 10**5)*(t/10**6)**19.9



energycno = []
for i in range(0,len(T)):
    a = ecno(rho[i],T[i])
    energycno.append(a)

fig = plt.figure()
plt.plot(rRstar,energycno,color='salmon')
plt.title("Energy generation rate (CNO) versus radius" )
#plt.xticks(np.arange(0, 10, 1))
plt.ylabel("Energy generation rate")
plt.xlabel('r/Rsun')
plt.xlim(0,1)
plt.grid()
#fig.savefig('3i2.png')
plt.show()


dLdr = []

for i in range(0,len(rho)):
    dL = 4*pi*(rRstar[i])**2 *rho[i]*(energycno[i] + energypp[i])
    dLdr.append(dL)

fig = plt.figure()
plt.plot(rRstar,dLdr,color='salmon')
plt.title("dL/dr versus radius" )
#plt.xticks(np.arange(0, 10, 1))
plt.ylabel("dL/dr")
plt.xlabel('r/Rsun')
plt.xlim(0,1)
plt.grid()
#fig.savefig('3i3.png')
plt.show()


luminosity = np.trapz(dLdr[:650],rRstar[:650])
print luminosity 
print luminosity * Rsun**3 / (3.828*10**26)
print T[0]
print rho[0]



