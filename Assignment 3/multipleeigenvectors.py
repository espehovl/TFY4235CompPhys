# -*- coding: utf-8 -*-
"""
Created on Fri May  1 16:58:53 2020

@author: Espen
"""
import matplotlib.pyplot as plt
import numpy as np
from utilities import normalized_function, normalizationFactor

def Reader(filename): #Read from file - each discretization has its own file
	file = open(filename,"r")
	raw = file.readlines()
	x_space = [float(x) for x in raw[0].split()]
	vectors = []
	for i in range(1,len(raw)):
		vectors.append([float(r) for r in raw[i].split()])
	file.close()
	return x_space, vectors

def analytical(x,n):
	return np.sqrt(2)*np.sin(n*np.pi*x)

fig = plt.figure()
#Calculate errors for each eigenvector
norm = normalizationFactor(Reader("results/vectors/vector10.tsv")[0],Reader("results/vectors/vector10.tsv")[1][0])
MSE = [] #Mean squared errors for the 10 different eigenvalues, for each discretization
for i in range(10,6111,100): #Go through all different discretizations
	discMSE = [] #The mean squared errors for the current discretization
	x,v = Reader("results/vectors/vector"+str(i)+".tsv")
	vectorIndex = 0
	for vector in v: #For each eigenvector
		analytic = [analytical(X,vectorIndex+1) for X in x]
		temp_errors = [] #The squared errors for the current vector
		for j in range(len(vector)): #For each value in the eigenvector
			temp_errors.append((norm*vector[j]-analytic[j])**2)
		discMSE.append(sum(temp_errors)/len(temp_errors))
		plt.scatter(i,discMSE[-1])
		vectorIndex += 1
	MSE.append(discMSE)

		

plt.legend([r"$\psi_{{{}}}$".format(i) for i in range(1,11)])
plt.xlabel("N")
plt.ylabel("MSE")
#plt.xscale("log")
#plt.yscale("log")
plt.show()