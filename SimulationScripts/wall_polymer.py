#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:37:49 2019

@author: gaurav
	 Greg 
"""

### BOX SIZE -100 100 -100 100 -100 100
#50 bead polymer

import sys
import numpy as np
import random
try:
    org_Type = str(sys.argv[1])
    num_Hinge = int(sys.argv[2])
    wall = 0
    seed = int(sys.argv[4])
except:
    print("ERROR: Unable to Parse Input Arguments\n - Setting Hinge Number to Zero")
    org_Type = "n"
    num_Hinge = 0
    wall = 0
    seed = 10

sub_Wall = 0
poly_length = 50
num_atoms = poly_length

if(wall == 1):
    num_atoms= num_atoms +20000
if(sub_Wall==1):
    num_atoms = num_atoms + 19602

file_object=open("initial.xyz","w+")
file_object.write("\n%d atoms \n%d bonds\n%d angles" %(num_atoms, poly_length-1, poly_length-2-num_Hinge))
file_object.write("\n\n2\t atom types\n1\t bond types\n1\t angle types")
file_object.write("\n\n-50 50 xlo xhi\n-50 50 ylo yhi\n-51 51 zlo zhi")
file_object.write("\n\nMasses \n\n1\t1.0 \n2\t1.0")
file_object.write("\n\nAtoms \n \n")



i=1
##Put down monomer beads
for x in np.arange(-50,poly_length-50,1):
    file_object.write("\t%d\t%d\t%d\t%0.2f\t%0.2f\t%0.2f\n" %(i,1,1,x,0,0.0))
    i+=1

if (wall ==1):
##Put down wall beads (Upper wall)
    for x in np.arange(-49.5,50,1):
        for y in np.arange(-49.5,50,1):
            file_object.write("\t%d\t%d\t%d\t%0.2f\t%0.2f\t%0.2f\n" %(i,2,2,x,y,50.0))
            i+=1
##Put down wall beads (Lower wall)
    for x in np.arange(-49.5,50,1):
        for y in np.arange(-49.5,50,1):
            file_object.write("\t%d\t%d\t%d\t%0.2f\t%0.2f\t%0.2f\n" %(i,2,2,x,y,-50.0))
            i+=1

if(sub_Wall==1):
    ##Put down wall beads (Upper wall)
    for x in np.arange(-49,50,1):
        for y in np.arange(-49,50,1):
            file_object.write("\t%d\t%d\t%d\t%0.2f\t%0.2f\t%0.2f\n" %(i,3,3,x,y,50.0))
            i+=1
    ##Put down wall beads (Lower wall)
    for x in np.arange(-49,50,1):
        for y in np.arange(-49,50,1):
            file_object.write("\t%d\t%d\t%d\t%0.2f\t%0.2f\t%0.2f\n" %(i,3,3,x,y,-50.0))
            i+=1

file_object.write("\nBonds \n \n")
for i in np.arange(1,poly_length):
    file_object.write("\t%d\t%d\t%d\t%d\n" %(i,1,i,i+1))
if (poly_length-2-num_Hinge) != 0 :
    file_object.write("\nAngles \n \n")
j = 0
bondFill = range(1,poly_length-1)

if(org_Type == "r" or org_Type=="u"):
    random.seed(seed+num_Hinge)
    bondFill = random.sample(bondFill,k=(poly_length-2-num_Hinge))
    spacing = (poly_length-2)/(num_Hinge)
    extraSpace = np.zeros(10)
    for k in range(int((np.rint((spacing-np.floor(spacing))*10)))):
        extraSpace[k] = 1
    np.random.shuffle(extraSpace)
    last = 0


nIncrease = 0
for i in np.arange(1,poly_length-1):
    if(org_Type=="n"):
        if( i > (poly_length-num_Hinge)/2-1 and (i<np.floor((poly_length-num_Hinge)/2) +num_Hinge)):
            j += 1
            continue
    elif( org_Type == "c"):
        if(i<=num_Hinge):
            j += 1
            continue
    elif (org_Type == "r"):
        if(i not in bondFill):
            continue
    elif ( org_Type=="u"):

        spacing = (poly_length-2)/(num_Hinge)
        if i-last ==(np.floor(spacing)+extraSpace[j % 10]):
            j = j+1
            last = i
            continue


    file_object.write("\t%d\t%d\t%d\t%d\t%d\n" %(i-j,1,i,i+1,i+2))

file_object.close()
