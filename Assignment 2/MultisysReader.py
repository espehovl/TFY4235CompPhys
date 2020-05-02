# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 20:37:44 2020

@author: Espen
"""
import pickle
print("Reading data from file...")

##Read data from simulations 
## This allows for reading large datafiles and saving their contents as Python lists. convenient

#file = open("C:/Users/Espen/OneDrive – NTNU/TFY4235 Numfys/Assignment 2/Results/Task 13/multisysPART1000particlesPotentialOn.txt","r")
file = open("C:/Users/Espen/OneDrive – NTNU/TFY4235 Numfys/Assignment 2/multisys2.txt","r")

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
print("Data loaded!")
print("Dumping pickles...")

#Write data to binary files - reading files takes a lot of time!
pickle.dump(T_mult,open("t.txt","wb"))
pickle.dump(X_mult,open("x.txt","wb"))

print("Pickles are pickled!")

del T_mult
del X_mult
del T_temp
del X_temp

