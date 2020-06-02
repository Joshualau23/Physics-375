# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 14:12:48 2019

@author: joshu
"""
import numpy as np
from math import log, pi
import matplotlib.pyplot as plt
import datetime
from astropy.io import ascii
from scipy.stats import linregress
from astropy import constants


###Q2a


data=ascii.read('hipparcos.txt')
parallax1 = data['col1']
parallax = parallax1*(0.001)
Vmag = data['col2']
Bmag = data['col3']
Imag = data['col4']

distance = []

for i in range(0,len(parallax)):
    d = 1.0 / parallax[i]
    distance.append(d)
    
Magv = []
for k in range(0,len(distance)):
    absmag = Vmag[k] - 5*log(distance[k] / 10.0,10)
    Magv.append(absmag)

BV = Bmag-Vmag

fig,ax=plt.subplots()
ax.scatter(BV,Magv)
plt.ylabel('Absolute Magnitude V')
plt.xlabel('B-V')
plt.gca().invert_yaxis()
plt.yticks(np.arange(-7, 18, 1.5))
plt.title('Absolute Magnitude V vs B-V ')
plt.show()

###Q2b


logtemperature = []

for i in range(0,len(BV)):
    temp = 9000.0 / ((BV[i]) + 0.93)
    a = log(temp,10)
    logtemperature.append(a)

logluminosity = []

for k in range(0,len(Magv)):
    lum = 10**(0.4*(4.83 - Magv[k]))
    a = log(lum,10)
    logluminosity.append(a)

fig,ax=plt.subplots()
plt.scatter(logtemperature,logluminosity)
plt.ylabel('log(Lv/Lsun)')
plt.xlabel('log(T)')
plt.title('log(Lv/Lsun) vs log(T) ')
plt.show()

###Q2c


steboltz =  5.670*(10**(-8))
def lumino(R,T):
    return 4*pi*steboltz*(R**2)*((T)**4)

Lsun = 3.828*(10**26)
rsun = 695700000
luminosity1 = []
luminosity02 = []
luminosity5 = []
temperature = []

for i in range(0,len(BV)):
    temp = 9000.0 / ((BV[i]) + 0.93)
    temperature.append(temp)


for j in range(0,len(temperature)):
    lum1 = log((lumino(1*rsun,temperature[j])) / Lsun ,10)
    lum02 = log((lumino(0.2*rsun,temperature[j])) / Lsun,10)
    lum5 = log((lumino(5*rsun,temperature[j])) / Lsun,10)
    luminosity1.append(lum1)
    luminosity02.append(lum02)
    luminosity5.append(lum5)
    
    
    
plt.scatter(logtemperature,luminosity1,color='blue')
plt.scatter(logtemperature,luminosity02,color='green')
plt.scatter(logtemperature,luminosity5,color='r')
plt.ylabel('log(L/L0)')
plt.xlabel('log t')
plt.title('log(Lv/Lsun) vs log(T) ')
plt.show()


logtemperature = []

for i in range(0,len(BV)):
    temp = 9000.0 / ((BV[i]) + 0.93)
    a = log(temp,10)
    logtemperature.append(a)

logluminosity = []

for k in range(0,len(Magv)):
    lum = 10**(0.4*(4.83 - Magv[k]))
    a = log(lum,10)
    logluminosity.append(a)

fig,ax=plt.subplots()
plt.scatter(logtemperature,logluminosity)
plt.plot(logtemperature,luminosity1,color='blue')
plt.plot(logtemperature,luminosity02,color='green')
plt.plot(logtemperature,luminosity5,color='r')
plt.ylabel('log(Lv/Lsun)')
plt.title('log(Lv/Lsun) vs log(T) ')
plt.xlabel('log(T)')
plt.show()


### q3a

data=ascii.read('W19_assignment1_orbit.dat')
orbitalphase =  data['col1']
radialvelocity1 = data['col2']
radialvelocity2 = data['col3']
spectro = data['col4']

time = []
for i in range(0,len(orbitalphase)):
    t = 50.0 * orbitalphase[i]
    time.append(t)
    

fig,ax=plt.subplots()
ax.scatter(time,radialvelocity1)
ax.scatter(time,radialvelocity2)
plt.ylabel('Radial Velocity (km/s)')
plt.xlabel('Time (Days)')
plt.title('Radial Velocities vs Time ')
plt.show()


###q3b

radialvelocity1round = np.around(radialvelocity1,decimals = 1)
radialvelocity2round = np.around(radialvelocity2,decimals = 1)
samevelocity = []

for i in range(0,len(radialvelocity1round)):
    if radialvelocity1round[i] == radialvelocity2round[i]:
        samevelocity.append(radialvelocity1round[i])
    else:
        continue


maxvr1 = max(radialvelocity1round) - samevelocity[0]
maxvr2 = max(radialvelocity2round) - samevelocity[0]

P = 4320000
def msin3i(n,m):
    return (P / (2*pi*G))*((1+ (n/m))**-1)*((n+m)**3)

M1 = msin3i(maxvr1,maxvr2) * 1000**3
M2 = msin3i(maxvr2,maxvr1)* 1000**3
print M1,M2

time2 = []
for i in range(0,len(time)):
    t2 = time[i] + 52.0
    time2.append(t2)
    
    


fig,ax=plt.subplots()
ax.scatter(time,spectro, color = 'gold')
ax.scatter(time2,spectro, color = 'gold')
plt.ylabel('Apparent Magnitude')
plt.xlabel('Time (Days)')
plt.gca().invert_yaxis()
plt.title('Apparent Magnitude vs Time ')
plt.show()

###q3c

m0 = min(spectro)

def luminoratio(m,m0):
    return 100**((m-m0) / 5.0)


loglml0ratio = []

for i in range(0,len(spectro)):
    a = log(luminoratio(spectro[i],m0),10)
    loglml0ratio.append(a)


fig,ax=plt.subplots()
ax.scatter(time,loglml0ratio)
plt.ylabel('Log(L/L0)')
plt.xlabel('Time (Days)')
plt.title('log(L/L0) vs Time ')
plt.show()


###q3d

mld = max(spectro)
msd = 1.507312


#for i in range(0,len(spectro)):
 #   if spectro[i] > 1.5 and spectro[i] < 1.517 :
  #      msd.append(spectro[i])
   # else:
    #    continue


#msd.sort()
#for i in range(0,len(msd) - 1):
 #   if round(msd[i],4) == round(msd[i+1],4):
  #      msd2.append(msd[i])
   # else:
    #    continue

#print msd2

#rh / rc
tempratio = (1.0 - luminoratio(mld,m0)) / (1.0 - luminoratio(msd,m0))
print tempratio



###q3d

ta = 0.473526
tb = 0.47952
tc = 0.52048

#ta1 = []
#tb1 = []

#for i in range(0,len(spectro)):
#    if spectro[i] == mld and orbitalphase[i] < 0.6 and orbitalphase[i] > 0.4:
#        tb1.append(orbitalphase[i])
        
#print tb1
        
#rs / rl
radiusratio = (tb - ta) / (tc - ta)
print radiusratio







