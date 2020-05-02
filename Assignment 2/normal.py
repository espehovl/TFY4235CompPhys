# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 15:38:45 2020

@author: Espen
"""

import matplotlib.pyplot as plt

nums = []
file = open("normaldistribution.txt", "r")
for line in file.readlines():
	nums.append(float(line.split()[0]))
file.close()



fig = plt.figure()
plt.hist(nums, bins=90)
plt.title(r'Random numbers from number generator, $\mu$=0.0, $\sigma$ = 1.0')

