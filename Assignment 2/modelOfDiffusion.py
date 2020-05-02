# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 17:52:47 2020

@author: Espen
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"]=(11,11)
plt.rcParams["font.size"]=50


def heatmap2d(arr , title, limits, Cmap = 'magma'):
    fig,ax = plt.subplots()
    plt.imshow(arr, cmap = Cmap,extent=limits, aspect = "auto")
    plt.gca().invert_yaxis()
    plt.clim(0,85)
    #plt.colorbar()
    plt.title(title)
    plt.xlabel(r"x [$\mu$m]")

    plt.ylabel("t [s]")
    plt.tight_layout()
    plt.show()

def n(x,t, r_i):
    N = 1001
    kBT = 26.0e-3*1.609e-19 #Joule
    etta = 1.0e-3 #Pascal seconds
    gamma_i = 6 *np.pi*etta*r_i # kg/s
    D = kBT/gamma_i # m^2/s
    return N/(np.sqrt(4*np.pi*D*t*1e12))*np.exp(-((x)**2)/(4*D*t*1e12))

r1 = 12e-9 #m
r2 = 3*r1

timespan = np.linspace(1e-5,7,101) #seconds
positionspan = np.linspace(-60,60,100) #micrometers
positionspan2 = np.linspace(-33,33,100) #micrometers

##Using the diffusion formula v
arr_particle1 = np.array([[n(x,t,r1) for x in positionspan] for t in timespan])
arr_particle2 = np.array([[n(x,t,r2) for x in positionspan2] for t in timespan])

heatmap2d(arr_particle1,"",[-60,60,7,0])
heatmap2d(arr_particle2,"",[-33,33,7,0])


