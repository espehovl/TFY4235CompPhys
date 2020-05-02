# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 23:48:50 2020

@author: Espen
"""
import matplotlib.pyplot as plt



print("Plot coming up!")


plt.rcParams["figure.figsize"]=(5,4)
##Fetch position data from a single particle
file = open("pos_data.txt","r")
X = []
times = []
x_temp = []
t_temp = []
deltaT=[]
Tau=[]
for line in file.readlines(): #Read the data file lines
    try:
        if ((not line.split()[0] == "x") or (not line[0]=="#") ):
            x_temp.append(float(line.strip().split()[0]))
            t_temp.append(float(line.strip().split()[1]))
    except:
        if line[0] =="#":
            X.append(x_temp)
            times.append(t_temp)
            x_temp = []
            t_temp = []
            deltaT.append(float(line[1:].split()[0]))
            Tau.append(float(line[1:].split()[1]))
        else: 
            None
file.close()

##Fetch metadata
file = open("metadata.txt","r")
temp = file.readlines()[1].split()

dt=float(temp[0]) #s
tau = float(temp[1]) #s
L = float(temp[2]) # um
steps = int(temp[3])
R = float(temp[4]) #nm
deltaU = float(temp[5]) #eV
kBT = float(temp[6]) #eV

file.close()

##Fetch potential profile
file = open("pot.txt","r")
pots = []
X_pots = []

for line in file.readlines():
    pots.append(float(line.split()[0]))
    X_pots.append(float(line.split()[1]))
    
#Normalize the potential values, so that it matches the timescale better
m = max(pots)
if max(pots) == 0: m = 1
pots_norm=[i*max(times[0])/m for i in pots]

file.close()

##Plot the single particle

fig = plt.figure()
for i in range(len(X)):
    plt.plot(X[i], times[i], label=r'$\delta T=$'+str(deltaT[i]) + r', $\tau=$'+str(Tau[i]))

##Plot potential profile over the position tracks
plt.plot(X_pots,pots_norm, label="Potential profile, \nnot to scale!", color = "C1", alpha = 0.1)
plt.fill_between(X_pots,pots_norm, color="C1", alpha = 0.4)

## Display horizontal lines according to the period
#i = 0
#while (i*tau < max(pots_norm)):
#    plt.hlines(i*tau+0.75*tau,min(X_pots), max(X_pots), alpha = 0.5, linestyle="dashed")
#    plt.hlines(i*tau+tau,min(X_pots), max(X_pots), alpha = 0.5, linestyle="dashed")
#    i+=1

##plt.suptitle(f"Random walk of particle in a sawtooth potential, {steps} steps")
#plt.title(r"$\delta$T="+str(dt)+r"s, $\tau=$"+str(tau)+r"s, L="+str(L)+f"µm, R="+str(R)+ r"nm, $\Delta$U="+str(deltaU)+r"eV, $k_B T =$"+str(kBT)+"eV", size=10)

plt.title(r"$\delta$t="+str(dt)+r" s, L="+str(L)+f" µm, $\Delta$U="+str(deltaU)+r" eV")
plt.ylabel(r'$t$ [s]')
plt.xlabel(r'$x$ [$\mu$m]')
#plt.legend()
fig.savefig("RandomWalk.png")
plt.show()


##Make a histogram of all x positions
"""
absPos = [abs(i) for i in X[0]]
fig = plt.figure()
plt.hist(absPos,bins=60)
plt.title("Histogram of positions")
plt.xlabel("|x| [µm]")
fig.savefig("PosHistogram.png")
plt.show()

"""


#%%
##Plot multiparticle data
print("Multisystem plots coming up!")

file = open("MetadataMultisys.txt","r")

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
#Get the average drift velocity

driftAvg = sum(driftVelsMult)/len(driftVelsMult)
print("Average drift velocity:", driftAvg)

## Plot multiple trajectories at once
##Read data
file = open("multisys2.txt","r")
"""
T_mult = [0]
X_mult = [[float(i)] for i in file.readline().split()[1:]]
for line in file.readlines():
    ll = line.split()    
    T_mult.append(float(ll[0]))
    for i in range(len(ll[1:])):
        X_mult[i].append(float(ll[i+1]))
"""
T_mult = []
X_mult = []

T_temp = []
X_temp = []
for line in file.readlines():
	ll = line.split()
	if ll[0] != "#":
		T_temp.append(float(ll[0]))
		X_temp.append(float(ll[1]))
	elif ll[0] == "#":
		T_mult.append(T_temp)
		X_mult.append(X_temp)
		T_temp = []
		X_temp = []


fig=plt.figure()
for ii in range(len(X_mult)):
    plt.plot(X_mult[ii],T_mult[ii])#, label=r"$\tau = $"+str(tauMult[ii])+" s")
    
plt.title(f"Positions of {len(X_mult)} simulations")
plt.xlabel("x [µm]")
plt.ylabel("t [s]")
plt.legend(["Small particle","Large particle"])

plt.plot(X_pots,pots_norm, label="Potential profile, \nnot to scale!", color = "C1", alpha = 0.1)
plt.fill_between(X_pots,pots_norm, color="C1", alpha = 0.4)

fig.savefig("multisystemPos.png")
plt.show()


#%%
##Locate and plot longest walker
"""
maxPos = max(endPosMult)
maxIndex = endPosMult.index(maxPos)
fig = plt.figure()

plt.plot(X_mult[maxIndex], T_mult)
plt.legend([r"$\tau$ = "+str(tauMult[maxIndex])])
plt.title("Trajectory of longest walk")
plt.xlabel("x [µm]")
plt.ylabel("t [s]")

fig.savefig("maxPos.png")

plt.show()
"""
##Plot the drift velocity against the time of ratchet being off
"""
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
#%%
z = tauMult[::degeneracy]
#%%

print(f"Max velocity is {max(driftVelAvg)} µm/s, which occurs at tau = {z[driftVelAvg.index(max(driftVelAvg))]}")

plt.plot(z,driftVelAvg)
plt.title("Drift velocity vs. time periodicity of ratchet")
plt.xlabel(r"$\tau$ [s]")
plt.ylabel(r"$\langle v_{drift} \rangle$ [µm/s]")
fig.savefig("driftvelocities.png")
plt.show()

"""

