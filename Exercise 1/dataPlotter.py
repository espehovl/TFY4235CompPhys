# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 16:19:18 2020

@author: Espen
"""

import matplotlib.pyplot as plt

file = open("baksnepperdata.txt","r")

species = []
ages = []

raw = file.readlines()
for elem in raw:
	vals = elem.split(",")[0].split()
	age = elem.split(",")[1].split()
	temp = []
	for V in vals:
		temp.append(int(V))
	species.append(temp)
	del vals
	del temp
	temp = []
	for A in age:
		temp.append(int(A))
	ages.append(temp)
	del age
	del temp

file.close()

deathsonly = []
for el in ages:
	row = []
	for i in range(len(el)):
		if el[i]==0: row.append(1)
		else: row.append(0)
	deathsonly.append(row)


plt.style.use("seaborn-muted")

fig = plt.figure()
plt.title("Evolution progress")
plt.set_cmap("RdGy")
plt.pcolormesh(species, vmin = 0, vmax = 99, snap = True)
plt.ylabel("Duration")
plt.xlabel("Species #")

fig2 = plt.figure()
plt.title("Mutation activity")
plt.set_cmap("Greys")
plt.pcolormesh(deathsonly, snap = True)
plt.ylabel("Duration")
plt.xlabel("Species #")