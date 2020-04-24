# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 17:19:19 2020

@author: Espen
"""

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})
from utilities import readEigenData, normalized_function, findAlphas, innerProduct, animateFunction


disc,vals,vecs,potential, v0 = readEigenData("results/eigenDataBarrier3000.tsv") #Discretization, eigenvalues, eigenvectors, potential profile

#%%
numGraphs = 3 #How many eigenfunctions to show in one plot?

#Exact solution, properly normalized
analyticSol = lambda n,x: np.sqrt(2)*np.sin(np.pi*n*x)

analyticalFunctions = []
selectedNumFunctions = []
numericalNormalized = []


##Create lists of wavefunctions
for n in range(min(numGraphs,len(vals))):
	analyticalFunctions.append([analyticSol(n+1,x) for x in disc])
	selectedNumFunctions.append([(el*vals[n]) for el in vecs[n]])
	numericalNormalized.append(normalized_function(disc,vecs[n]))


## Create lists of wavefunctions squared
analyticalSquared = []
numericalSquared = []
errorsOfSquaredFuncs = []
for n in range(min(numGraphs,len(vals))):
	analyticalSquared.append([i**2 for i in analyticalFunctions[n]])
	numericalSquared.append([i**2 for i in numericalNormalized[n]]) 	
	errorsOfSquaredFuncs.append([analyticalSquared[n][i]-numericalSquared[n][i] for i in range(len(analyticalSquared[n]))])

potential_norm = [int(bool(v0))*i*max(numericalSquared[0]) for i in potential]
fig, (ax1) = plt.subplots(1,figsize = (6,4))

plt.tight_layout()
for ax in fig.get_axes():
    ax.label_outer()
#Plot numerical functions
#ax1.plot(disc,potential_norm,  color = "k")
ax1.fill_between(disc,potential_norm, color = "k", alpha = 0.5)
for i in range(min(numGraphs,len(vals))):
	ax1.plot(disc,numericalNormalized[i], label = r"$\psi_{{{}}}$".format(i+1))

#ax1.set_title(r"Numerical solution, normalized, $\nu_0 = $"+f"{v0}")
ax1.set_xlabel("x/L")
ax1.set_ylabel(r"$|\psi(x')|^2$")
ax1.legend()

fig.savefig("plot.png")
plt.show()
"""	
#Plot analytical functions
for i in range(min(numGraphs,len(vals))):
	ax2.plot(disc,analyticalSquared[i], label = r"$\psi_{{{}}}$".format(i+1))

ax2.set_title("Analytical solution")
ax2.set_ylabel(r"$|\psi(x')|^2$")
ax2.legend()

#Plot error in squared functions
for i in range(min(numGraphs,len(vals))):
	ax3.plot(disc,errorsOfSquaredFuncs[i], label = r"$\epsilon_{{{}}}$".format(i+1))
ax3.set_ylabel("Error")
ax3.set_xlabel("x/L")
ax3.legend()
plt.show()
"""
### Plot difference (ratio) between analytical and numerical solution
#fig = plt.figure(figsize=(6,4))
#for i in range(min(numGraphs, len(vals))):
#	plt.plot(disc,[abs(analyticalFunctions[i][idx])/abs(numericalNormalized[i][idx]) for idx in range(len(analyticalFunctions[i]))])
#plt.title("Ratio of solutions")
#plt.show()

#%%
#2.5 Plot computed eigenvalues against analytical solution eigenvals
"""
analyticalEigVals = [(np.pi * i)**2/(np.pi)**2 for i in range(1,len(vals)+1)]
numericalEigVals = [vals[i]/vals[0] for i in range(len(vals))]

E_0 = lambda m,h_bar,L: 2*m/(h_bar**2*L**2) #Relation found in 2.3


n_space =  [i for i in range(1,len(vecs[0])+1)]


fig = plt.figure(figsize = (6,4))

plt.plot(n_space[:len(analyticalEigVals)], analyticalEigVals, label = r"Analytical solution, $\lambda _n = (\pi n)^2$")
plt.plot(n_space[:len(numericalEigVals)], numericalEigVals, label = r"Numerical solution, $\lambda_n = \dfrac{2m}{\hbar^2L^2} E_n$")
plt.ylabel(r"$\lambda_n (E_n) / E_0$")
plt.xlabel(r"$n$")
plt.title("Eigenvalues as a function of n")
plt.legend()
plt.show()
"""

#%%

#2.8 calculate alpha_n

alphaN = findAlphas(numericalNormalized,disc)
print(alphaN)



#%%

#2.9 orthogonality check
print(innerProduct(numericalNormalized[1],numericalNormalized[1],disc))


#%%

##2.10 check initial condition
plt.rcParams.update({'font.size': 14})
#The real deal function, returns the time evolution of the entire system
def Psi(t_values, alphaN, lambdaN, psiN, initialPsi):
	"""
	t_values: t-space (list)
	
	The following lists must have equal lengths vvv
	alphaN: values for alpha at each n (list)
	lambdaN: eigenvalues for each n (list)
	psiN: eigenvectors for each n (2D-list)
	"""
	resultMat = np.zeros((len(t_values),len(psiN[0])),dtype=np.complex128) #The results at different times. Each row is the superposition of all wavefunctions psiN

	#Apply initialcondition
#	for i in range(resultMat.shape[1]): #For each element in the first row...
#		resultMat[0,i] = initialPsi[i]
	for tt_idx in range(0,len(t_values)): 	#For each time t...
		for n in range(len(psiN)): 	 	 	#for each eigenvector...
			for i in range(len(psiN[n])): 	#For each value of that eigenvector...
				resultMat[tt_idx,i] += psiN[n][i]*alphaN[n]*np.exp(-1j*lambdaN[n]*t_values[tt_idx])
	return resultMat
	
t_end = np.pi/(vals[1]-vals[0]) #End time

t = np.linspace(0,1,100) #Time space

#eigenFuncs = [[np.sqrt(2)*np.sin(n*np.pi*x) for x in disc] for n in range(1,6)]
eigenFuncs = selectedNumFunctions
Psi_0 = [1/np.sqrt(2)*(eigenFuncs[0][i] + eigenFuncs[1][i]) for i in range(len(eigenFuncs[0]))]


#Psi_0 = np.zeros(len(Psi_0))
#Psi_0[len(Psi_0)//2] = 1
#Psi_0 = normalized_function(disc,Psi_0) #Dirac delta?

alpha_0 = [innerProduct(np.conj(V),Psi_0,disc) for V in eigenFuncs]
lambda_0 = vals[:len(eigenFuncs)]

basicTest = Psi(np.multiply(t,np.pi/(vals[1]-vals[0])), alpha_0, lambda_0, eigenFuncs, Psi_0) #Make a 2D array of superpositions of psi's over time

## Plot the behemoth - plot the superposition of wavefunctions as time progresses

fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize = (12,7))
for ax in fig.get_axes():
    ax.label_outer()

#Plot real parts
reals = [[n.real for n in bT] for bT in basicTest ]

C1 = ax1.contourf(t,disc,np.transpose(reals)) #Transpose is to let time be plotted along x-axis and position along y-axis
fig.colorbar(C1, ax = ax1)
ax1.set_title(r"$\Re[\Psi(x',t')]$")
ax1.set_ylabel("x'")

#Plot imaginary parts
imags = [[n.imag for n in bT] for bT in basicTest ]
C2 = ax2.contourf(t,disc,np.transpose(imags))
fig.colorbar(C2, ax = ax2)
ax2.set_title(r"$\Im[\Psi(x',t')]$")
ax2.set_ylabel("x'")

#Plot squared wavefunctions
squareds = [[np.square(np.abs(n)) for n in bT] for bT in basicTest ]
C3 = ax3.contourf(t,disc,np.transpose(squareds))
fig.colorbar(C3, ax=ax3)
ax3.set_title(r"$|\Psi(x',t')|^2$")
ax3.set_xlabel(r"$t'\ \dfrac{\pi}{\lambda_2-\lambda_1}$")
ax3.set_ylabel("x'")

plt.show() 

#%% Animate the tunneling (mostly for fun)

#plt.ioff()
#animateFunction(squareds, disc, potential)



#%% Task 3.4


def f(_lambda):
	k = np.sqrt(_lambda)
	kappa = np.sqrt(v0-_lambda)
	return np.exp(kappa/3)*np.square(kappa*np.sin(k/3)+k*np.cos(k/3)) - np.exp(-kappa/3)*np.square(kappa*np.sin(k/3)-k*np.cos(k/3))

def df(_lambda): #This looks horrible, but it is correct! (or so it seems!)
	k = np.sqrt(_lambda)
	kappa = np.sqrt(v0-_lambda)
	
	return (-1)/6*(k*np.cos(1/3*k)+kappa*np.sin(1/3*k))**2 * np.exp(1/3*kappa)/kappa - 1/6*(k*np.cos(1/3*k)+kappa*np.sin(1/3*k))**2 * np.exp((-1)/3*kappa)/kappa + 2*(k*np.cos(1/3*k)+kappa*np.sin(1/3*k))*((-1)/2*np.sin(1/3*k)/kappa+1/2*np.cos(1/3*k)/k+1/6*kappa*np.cos(1/3*k)/k-1/6*np.sin(1/3*k))*np.exp(1/3*kappa)-2*(k*np.cos(1/3*k)+kappa*np.sin(1/3*k))*((-1)/2 *np.sin(1/3*k)/kappa+1/2*np.cos(1/3*k)/k+1/6*kappa*np.cos(1/3*k)/k-1/6*np.sin(1/3*k))*np.exp((-1)/3*kappa)

def findRoots(guesses, tolerance = 1e-5):
	"""
	guesses: a list of guesses of the roots we are looking into finding
	tolerance: a tolerance of how precise we want the roots
	"""
	roots = []
	iterations = [] #Mostly for debugging
	
	#Newton's method. Pretty basic, but works for finding separated roots. 
	#Not sure if I will spend more time on this
	for g in guesses: #For all guesses...
		if g > v0: #We are outside the scope of interest
			break
		xOld = g #Initial guess
		counter = 0
		while True: #Perform iterations
			counter += 1
			xNew = xOld - f(xOld)/df(xOld)
			if abs(xNew - xOld) < tolerance: #Solution is satisfactory
				break
			xOld = xNew
		iterations.append(counter)
		roots.append(xNew)
	print(roots, iterations, sep="\n")
	return roots

def bruteRoots(guess):
	#Brute force - iterate through function values, find where the sign changes.
	#Assumes/requires very fine discretization!
	roots = []
	dX = 1 #Delta, to avoid searching too large areas
	
	for g in guess: #For each guess...
		if g > v0: 
			break #We are outside the interesting area
		xSpace = np.linspace(g-dX, g+dX, 101) #Create a fine discretization xspace around the guess
		fSpace = [f(i) for i in xSpace] #Find the values...
		for i in range(len(fSpace)-1):
			if fSpace[i]*fSpace[i+1] < 0:
				roots.append(xSpace[i])
	print(roots)
	return roots

lambdaSpace = np.linspace(0,int(v0),int(v0*10)+1)
fVals = [f(L) for L in lambdaSpace]
dfVals = [df(L) for L in lambdaSpace]
fig = plt.figure(figsize=(6,4))
plt.plot(lambdaSpace,fVals)
plt.plot(lambdaSpace,dfVals)

plt.vlines(bruteRoots(vals), 0, max(fVals), linestyle = "dashed")
plt.hlines(0,0,v0)

plt.ylabel(r"$f(\lambda)$")
plt.xlabel(r"$\lambda$")
plt.show()



