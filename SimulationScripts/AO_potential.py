#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 01:34:01 2022

@author: gaurav
"""
import numpy as np
import matplotlib.pyplot as plt
import os

for os_p in [0.80, 0.90]:
    os.mkdir("p_"+str(os_p))
    file_object=open("p_"+str(os_p)+"/LJ_AO.table","w+")
    file_object.write("\nLJ_AO_T \nN 1901\tR\t0.10\t2.0\n\n")
    i=0
    pot,force,dist=[],[],[]
    for r in np.arange(0.10,2.002,0.001):
        dist.append(r)
        if(r<1):
            pot.append(4*1/(r**12)-4*1/(1**12)-np.pi*os_p*5/12)
            force.append(48/(r**13))
        if(r>1):
            pot.append(-np.pi*os_p/12*((2-r)*(2-r)*(4+r)))
            force.append(-np.pi*os_p/4*(4-r*r))
        i+=1
    print(r)
    
    print(i)
        
    for i in range(1,1902):
        file_object.write("%d\t %0.4f\t %0.4f\t %0.4f\n"%(i,dist[i-1],pot[i-1],force[i-1])) 
    
    file_object.close()
    
