# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 13:17:31 2019

@author: joshu
"""

import numpy as np
from numpy import arange,amax,amin,array
from math import log, pi,e
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.stats import linregress
from scipy.constants import speed_of_light,Planck,Boltzmann


###1a


def blackbodywave(lamda,temp):
    return ((2*Planck*speed_of_light**2) / lamda**5)*((e**(Planck*speed_of_light / (lamda*Boltzmann*temp)) - 1.0)**-1)

waverange = arange((100*(10**-9)),(2*(10**-6)),(5*(10**-10)))

BBwavelength = []

for i in range(0,len(waverange)):
    l = waverange[i]
    bb = blackbodywave(l,5500)
    BBwavelength.append(bb)
    
plt.scatter(waverange,BBwavelength)
plt.xlim((100*(10**-9)),(2.1*(10**-6)))
plt.xlabel('Wavelength')
plt.ylabel('Plank Function')
plt.title('Plank Function vs Wavelength')
plt.show()

print amax(BBwavelength)

peakwavelength = 0
for i in range(0,len(waverange)):
    if BBwavelength[i] == 20612853101415.254:
        peakwavelength = waverange[i]
    else:
        continue

print peakwavelength




###1b

def blackbodyfreq(mu,temp):
    return ((2*Planck*mu**3 / speed_of_light**2))*((e**(Planck*mu / (Boltzmann*temp)) - 1.0)**-1)

freqrange = arange((0.15*10**14),(30*10**14),6.5**14)

BBfrequency = []

for i in range(0,len(freqrange)):
    mu = freqrange[i]
    bb = blackbodyfreq(mu,5500)
    BBfrequency.append(bb)
    
plt.scatter(freqrange,BBfrequency)
plt.xlim((0.05*10**14),(31*10**14))
plt.ylim(0,(50*10**-9))
plt.xlabel('Frequency')
plt.ylabel('Plank Function')
plt.title('Plank Function vs Frequency')
plt.show()

print amax(BBfrequency)

peakfrequency = 0
for i in range(0,len(freqrange)):
    if BBfrequency[i] == 3.154577180833425*10**-8:
        peakfrequency = freqrange[i]
    else:
        continue

BBfrequency = array(BBfrequency)    
item_index = np.where(BBfrequency==3.154577180833425*10**-8)
print freqrange[item_index]





###1d

freqrange = arange((0.15*10**14),(30*10**14),6.5**14)

BBfreqtimefreq = []

for i in range(0,len(freqrange)):
    mu = freqrange[i]
    bb = mu*blackbodyfreq(mu,5500)
    BBfreqtimefreq.append(bb)


#plt.scatter(freqrange,BBfreqtimefreq)
#plt.xlim((0.05*10**14),(31*10**14))
#plt.ylim(0,(15*10**6))
#plt.show()

print amax(BBfreqtimefreq)


waverange = arange((100*(10**-9)),(4*(10**-6)),(15*(10**-10)))

logwavelength = []
logBBwavetimeswave5500 = []
logBBwavetimeswave3000 = []
logBBwavetimeswave30000 = []
logfrequency = []
logBBfreqtimefreq = []

for i in range(0,len(waverange)):
    l = waverange[i]
    logwav = log(waverange[i],10)
    logwavelength.append(logwav)
    
    bb5500 = log(l*blackbodywave(l,5500),10)
    logBBwavetimeswave5500.append(bb5500)
    bb3000= log(l*blackbodywave(l,3000),10)
    logBBwavetimeswave3000.append(bb3000)
    bb30000 = log(l*blackbodywave(l,30000),10)
    logBBwavetimeswave30000.append(bb30000)
    

    
plt.scatter(logwavelength,logBBwavetimeswave5500, color = 'aqua')
plt.scatter(logwavelength,logBBwavetimeswave3000,color = 'Gold')
plt.scatter(logwavelength,logBBwavetimeswave30000,color = 'salmon')
plt.xlabel('log(wavelength)')
plt.ylabel('log(Plank Function)')
plt.title('log(Plank Function) vs log(wavelength)')
plt.show()

viswaverange = arange((400*(10**-9)),(800*(10**-9)),(1*(10**-9)))

logviswavelength = []
logvisBBwavetimeswave5500 = []
logvisBBwavetimeswave3000 = []
logvisBBwavetimeswave30000 = []

for i in range(0,len(viswaverange)):
    l = viswaverange[i]
    logwav = log(viswaverange[i],10)
    logviswavelength.append(logwav)
    
    bb5500 = log(l*blackbodywave(l,5500),10)
    logvisBBwavetimeswave5500.append(bb5500)
    bb3000= log(l*blackbodywave(l,3000),10)
    logvisBBwavetimeswave3000.append(bb3000)
    bb30000 = log(l*blackbodywave(l,30000),10)
    logvisBBwavetimeswave30000.append(bb30000)

plt.scatter(logviswavelength,logvisBBwavetimeswave5500, color = 'aqua')
plt.scatter(logviswavelength,logvisBBwavetimeswave3000,color = 'Gold')
plt.scatter(logviswavelength,logvisBBwavetimeswave30000,color = 'salmon')
plt.xlabel('log(wavelength)')
plt.ylabel('log(Plank Function)')
plt.title('log(Plank Function) vs log(wavelength)')
plt.show()

vary30k = amax(logvisBBwavetimeswave30000) - amin(logvisBBwavetimeswave30000)
vary5k = amax(logvisBBwavetimeswave5500) - amin(logvisBBwavetimeswave5500)
vary3k = amax(logvisBBwavetimeswave3000) - amin(logvisBBwavetimeswave3000)

print vary30k,vary5k,vary3k


###3a

Rsun = 695500
Te = 10000
p = 10**-6 * 1000**3
k = 3.0 / 1000**2
def eddington(s):
    return ( (3.0/4.0) * Te**4 *( (p*k*s) + (2.0/3.0)))**(1.0/4.0)
def tau(s):
    return p*k*s

height = arange(0,1000,0.5)

temperature = []
ttau = []
for i in range(0,len(height)):
    x = eddington(height[i])
    t = tau(height[i])
    temperature.append(x)
    ttau.append(t)
    
plt.scatter(height,temperature)
plt.ylabel('Temperature')
plt.xlabel('Depth')
plt.title('Temperature vs Depth ')
plt.show()

#plt.scatter(height,ttau)
#plt.show()

ttau = array(ttau)
ttwooverthree = np.where(ttau==0.666 )
#print height[ttwooverthree]

temperature = array(temperature)

T10000 = np.where(temperature==9998.749765556617)
#print T10000
#print height[T10000]
#print height[T10000] / Rsun
#print (6.0/12.0)**(1.0/4.0)


####b



Kb = 8.6173303e-05
me = 9.10938356e-31
mp = 6.64465723e-27
kb = 1.38064852e-23
h = 6.62607004e-34
n = 10e19
rho = 10e-6


def saha(t):
    a = (mp / rho)*((2*pi*me*kb*t) / h**2)**(3.0/2.0) * exp(-13.6 / (Kb*t))
    return (a / 2.0)*(sqrt(1 + 4.0/a) - 1)

temp = arange(1,20000,50)
fraction = []

for i in range(0,len(height)):
   f = saha(temperature[i])
   fraction.append(f)

plt.scatter(temperature,fraction)
plt.ylim(0,1)
plt.show()


fs = []
for i in range(0,len(temperature)):
    f = saha(temperature[i])
    fs.append(f)

fs = array(fs)
plt.scatter(height,fs)
plt.ylim(0,1)
plt.ylabel('Fraction')
plt.xlabel('Depth')
plt.title('Fraction vs Depth ')
plt.show()

mostf2 = np.where(fs==0.9900108886802386)
#print height[mostf2]


###c

n = arange(3,103,1)

ang = 911.6e-10
def wavelength(n):
    return ang / ( 1 / 4.0 - 1.0/n**2) 

#print wavelength(3)


p = 10**-6 
k = 3.0 
kbal = 3.5e5  + k


height = arange(0,10,0.01)

def optdepth(s,kap):
    return p*kap*s


opticaldepthbalmer = []
opticaldepthreg = []
for i in range(0,len(height)):
    optbal = optdepth(height[i],kbal)
    optreg = optdepth(height[i],k)
    opticaldepthbalmer.append(optbal)
    opticaldepthreg.append(optreg)

#plt.scatter(height,opticaldepthbalmer)
plt.scatter(height,opticaldepthreg)
plt.ylim(0,1)
#plt.xlim(0,223)
plt.ylabel('Optical Depth')
plt.xlabel('Depth')
#plt.title('Optical Depth vs Depth ')
plt.show()

opticaldepthbalmer = array(opticaldepthbalmer)
tbalmer = np.where(opticaldepthbalmer==0.6650057)
#print height[tbalmer]

#print kbal / k

