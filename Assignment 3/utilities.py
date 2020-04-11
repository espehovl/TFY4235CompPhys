# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 17:04:20 2020

@author: Espen
"""

from scipy.integrate import simps
from scipy import sqrt
import numpy as np
import matplotlib.pyplot as plt

def readEigenData(filename):
	#Reads TSV (tab separated values) files with data of the form vv
	# Discretization
	# Potential v0 of barrier
	# Profile of potential barrier
	# Eigenvalue     Eigenvector
	discretization = []
	eigenVals = []
	eigenVecs = []
	potential = []
	v0 = 0
	file = open(filename,"r")
	
	counter = 0
	for line in file.readlines():
		counter +=1
		if line[0][0]=="#": 
			None
		elif counter == 2: #Fetch x_space
			discretization = [float(elem) for elem in line.split()]
		elif counter == 4: #Fetch potential barrier value v0
			v0 = float(line)
		elif counter == 6: #Fetch potential profile
			potential = [float(elem) for elem in line.split()]
		else: #Fetch eigenvalues and eigenvectors
			ll = line.split()
			eigenVals.append(float(ll[0]))
			eigenVecs.append([float(x) for x in ll[1:]])
	file.close()
	return discretization, eigenVals, eigenVecs, potential, v0

def average(_list):
	return sum(_list)/len(_list)

def normalized_function(x_values, y_values):
	integral = simps([np.square(np.abs(i)) for i in y_values],x_values)
	print("integral(absolute value) = ",np.abs(integral))
	return [i/np.sqrt(integral) for i in y_values]

def normalizationFactor(x_values,y_values):
	integral = simps([np.square(np.abs(i)) for i in y_values],x_values)
	print("integral(absolute value) = ",np.abs(integral))
	return 1/np.sqrt(integral)

##Find the alphas as according to task 2.8 and eqn. 2.14
##Eigenvectors are normalized
def findAlphas(eigenvectors, x_values):
	alphas = []
	for i in range(len(eigenvectors)):
		conjugateVec = np.conj(eigenvectors[i])
		integrand = np.multiply(conjugateVec,eigenvectors[i])
		alphas.append(simps(integrand, x_values))
	return alphas

def computeAlphaN(psiConjugate, Psi_0, x_values):
	return simps(np.multiply(psiConjugate, Psi_0),x_values)

def innerProduct(vec1, vec2, x_values):
	return simps(np.multiply(vec1,vec2), x_values)


def animateFunction(dataArray, x_space, potentialProfile, filename = "results/Animation/ani"):
	"""
	dataArray: a 2D-array of wavefunction superpositions as rows at different times, evolving over time
	x_space: the x_coordinates used
	potentialProfile: 1D-array of potentialProfile, for visual pleasure
	"""
	from PIL import Image
	import glob
	
	potentialProfile = np.multiply(potentialProfile,max(dataArray[1])+5)
	counter = 0
	plt.ioff()
	for row in dataArray:
		fig = plt.figure(figsize=(6,4))
		plt.fill_between(x_space,potentialProfile, color = "k", alpha = 0.5)
		plt.plot(x_space, row)
		plt.xlabel("x/L")
		plt.ylabel(r"$|\Psi(x',t')|^2$")
		fig.savefig(filename+f"{'{:04d}'.format(counter)}.png")
		counter += 1
		plt.close()
	
	# Create the frames
	frames = []
	imgs = glob.glob("results/Animation/*.png")
	for i in imgs:
		new_frame = Image.open(i)
		frames.append(new_frame)
		# Save into a GIF file that loops 	forever
		frames[0].save(filename+'.gif', format='GIF', append_images=frames[1:], save_all=True, duration=100, loop=1)

	print("GIF conversion complete!")