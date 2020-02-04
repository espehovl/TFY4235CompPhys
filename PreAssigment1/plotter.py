# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:13:01 2020

@author: Espen
"""

import matplotlib.pyplot as plt
import numpy as np

def fileReader():
	file = open("plotdata.txt","r")
	data = []
	for line in file.readlines():
		data.append((int(line.split(",")[0].strip()),int(line.split(",")[1].strip())))
	
	file.close()
	return data

data = fileReader()

xdata=[np.log(x[0]) for x in data]
ydata = [np.log(y[1]/1e6*100) for y in data]

coeffs = np.polyfit(xdata,ydata,deg = 1)
print("alpha = ",-coeffs[0])

plt.plot(xdata,ydata)


plt.xlabel("log #steps")
plt.ylabel("log #counts")
