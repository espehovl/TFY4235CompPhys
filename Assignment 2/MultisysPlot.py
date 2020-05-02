# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:48:43 2020

@author: Espen
"""

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"]=(7,4)


file = open("C:/Users/Espen/OneDrive – NTNU/TFY4235 Numfys/Assignment 2/Results/Task 9/MetadataMultisys0.01-5.0-100avgs.txt","r")

dtMult = []
tauMult = []
deltaUMult = []
kBTMult = []
endPosMult = []
driftVelsMult = []
degeneracy = int(file.readline().split()[-1])

for line in file.readlines():
    ll = line.split()
    dtMult.append(float(ll[1]))
    tauMult.append(float(ll[2]))
    deltaUMult.append(float(ll[6]))
    kBTMult.append(float(ll[7]))
    endPosMult.append(float(ll[8]))
    driftVelsMult.append(float(ll[9]))

file.close()

avgDriftVel = sum(driftVelsMult)/len(driftVelsMult)

print(avgDriftVel, "um/s")
##For plotting the average drift velocities as a function of tau

##Average over the degeneracy
fig = plt.figure()
driftVelAvg = []
#degeneracy = 100 #How many trials per choice of tau? 
tempSum = driftVelsMult[0]
for i in range(1, len(tauMult)+1):
    if i%degeneracy != 0:
        tempSum += driftVelsMult[i]
    elif i%degeneracy ==0:
        driftVelAvg.append((tempSum/degeneracy))
        if i < len(tauMult):
            tempSum = driftVelsMult[i]

z = tauMult[::degeneracy]


print(f"Max velocity is {max(driftVelAvg)} µm/s, which occurs at tau = {z[driftVelAvg.index(max(driftVelAvg))]}")

plt.plot(z,driftVelAvg)
plt.title("Drift velocity vs. time periodicity of ratchet")
plt.xlabel(r"$\tau$ [s]")
plt.ylabel(r"$\langle v_{drift} \rangle$ [µm/s]")

plt.show()
