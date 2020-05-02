# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 00:00:29 2020

@author: Espen
"""

import numpy as np
import matplotlib.pyplot as plt

print("Boltzmann snacks comin' up!")

def boltzmann(U, dU, kBT):
	return (np.exp(-U/kBT))/(kBT*(1-np.exp(-dU/kBT))) 



kBT = 26e-3 #eV
deltaU = 0.1*kBT #eV



Us = np.linspace(0,0.2)

boltz01 = [boltzmann(i, 0.1*kBT, kBT) for i in Us]
boltz10 = [boltzmann(i, 10 *kBT, kBT) for i in Us]

maxx = max(max(boltz01),max(boltz10))

b01norm = [i/maxx for i in boltz01]
b10norm = [i/maxx for i in boltz10]

##Plot the Boltzmann distribution

#plt.plot(Us, b01norm, label = r"$\Delta U = $"+str(0.1*kBT))
##plt.plot(Us, b10norm, label = r"$\Delta U = $"+str(10*kBT))
#plt.title("Boltzmann distribution")
#plt.legend()	
#plt.xlabel("U [eV]")
#plt.ylabel("p(U)")
#plt.show()


## Fetch end positions and plot histogram
#file = open("endpositions.txt","r")
#ends = []
#for line in file.readlines():
#	ends.append(abs(float(line)))
#
#fig = plt.figure()
#plt.hist(ends, bins = 30)
#plt.title("End positions, in absolute value")
#plt.xlabel("x [µm]")
#plt.show()
#fig.savefig("Endpositions.png")


##Fetch end position potentials and plot histogram
#file = open("potentialtracks.txt","r")
#Pots = []
#for line in file.readlines():
#	Pots.append(float(line))
#
#
#fig = plt.figure()
#plt.hist(Pots, bins = 120)
#plt.title("Potentials statistics")
#plt.xlabel("U [eV]")
#fig.savefig("PotentialStats.png")
#plt.show()