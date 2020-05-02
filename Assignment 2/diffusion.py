# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 19:51:22 2020

@author: Espen
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle

##Load binary files containing the lists of data (to keep from reading them over and over again)
print("Loadin' pickles...")
T_mult = pickle.load(open("t.txt","rb"))
X_mult = pickle.load(open("x.txt","rb"))
print("Pickles loaded!")


#Prepare for histogramming
print("Preparing for histogramming...")
X = []
T = []
for i in range(1,len(X_mult)): #Drop the first element, as it is from a different simulation
	for j in range(1,len(X_mult[1])):
		X.append(X_mult[i][j])
		T.append(T_mult[i][j])
		
print("Histogram list ready!")
#%%

plt.rcParams["figure.figsize"]=(10.5,11)
plt.rcParams["font.size"]=50


##Simulated for actual particle v

#Plot the particle densities over time
fig = plt.figure()
plt.hist2d(X,T, bins = 100, cmap = "magma")
plt.xlabel(r"x [$\mu$m]")
plt.xticks([0,40,80,120])
plt.ylabel("t [s]")
plt.clim(0,60000)
plt.tight_layout()
fig.savefig("particlehistogram.png")
#plt.colorbar()
plt.show()




