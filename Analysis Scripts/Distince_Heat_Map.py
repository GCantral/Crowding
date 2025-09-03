
import numpy as np
import matplotlib.pyplot as plt
import os,glob
import matplotlib as mpl


def FindHinges(poly_size):
    xyzFile = open("../initial.xyz", "r")
    line = ""
    while line !="Angles \n":
        line = xyzFile.readline()
        if not line:
            return range(2,poly_size)
    line = xyzFile.readline()
    line = xyzFile.readline()
    hinges = []
    for i in range(2,poly_size):
        if line and line.split()[3] == str(i):
            line = xyzFile.readline()
            continue
        hinges.append(i)
    return hinges


def CreateDIFile():
    filename = []
    for file in glob.glob("dump.*.txt"):
        filename.append(file)
    file= sorted(filename)
    total_steps=len(file)

    polymer_write = open("../DistInfo.txt", "w+")
    a = np.array(np.genfromtxt(file[0], skip_header=9))
    polysize  = len(a)
    hinges  = FindHinges(polysize)
    for i in hinges:
        polymer_write.write(str(i)+"\t")
    polymer_write.write("\n")
    dis_Total = np.zeros((polysize,polysize))
    k = 0
    for j in range(1000,total_steps,2):

        file_name=file[j]
        if(k==0):
            dis = np.zeros(polysize)
        a = np.array(np.genfromtxt(file_name, skip_header=9))
        a = np.array(sorted(np.array(a), key=lambda a_entry: a_entry[0]))        ### sort them by their monomer number
        for i in range(polysize):
            for z in range(i,polysize):
                dist = 0
                for k  in range(3):
                    dist = (pow(a[i,5+k]-a[z,5+k],2))+dist
                dis_Total[i,z] = dis_Total[i,z]+np.sqrt(dist)
                dis_Total[z,i] = dis_Total[i,z]


    dis_Total = 2*dis_Total/(total_steps-1000)
    for i in range(len(dis_Total)):
        line_str = str(dis_Total[i,0])
        for j in range(1,polysize):
            line_str = line_str+"\t"+str(dis_Total[i,j])
        polymer_write.write("%s\n"%(line_str))

    polymer_write.close()


def createInfo(parentFolder):
    os.chdir(os.getcwd()+"\\"+parentFolder)
    p_dirs = next(os.walk('.'))[1]
    for folder in p_dirs:
        os.chdir(os.getcwd()+"\\"+folder)
        i_dirs = next(os.walk('.'))[1]
        for i_folder in i_dirs:
            os.chdir(os.getcwd()+"\\"+i_folder+"\\coord_dump")
            CreateDIFile()
            print(os.getcwd())
            os.chdir("..")
            print(os.getcwd())
            os.chdir("..")
            print(os.getcwd())
        os.chdir("..")
    os.chdir("..")




def createDistinceHeatmap(parentFolder,Pic,hingeNum):
    os.chdir(os.getcwd()+"\\"+parentFolder+"\\p_0"+Pic+"\\K_50_"+hingeNum)
    readInfo = open("DistInfo.txt", "r")
    hinges = np.array(readInfo.readline().split())
    hinges = hinges.astype(int)
    readInfo.close()
    distData = np.array(np.genfromtxt("DistInfo.txt", skip_header=1))
    plt.rc('font',size = 15)


    for i in range (len(distData)):
        for j in range(i):
            distData[i][j]= distData[j][i]
            """
            if(distData[j][i] <3):
                distData[i][j] = 1
            else:
                distData[i][j]= 15"""


    plt.imshow(distData,cmap=mpl.colormaps["plasma_r"],vmin=0,vmax=15)
    ticks = ["0"]
    print(len(distData))
    if(int(hingeNum)==int(len(distData))-2):
        for i in range(int(hingeNum)-1):
            ticks.append("")
    else:
        for i in range(int(hingeNum)):
            ticks.append("")
    print(hinges)
    tickloc = [0]+list(hinges-1)+[len(distData)-1]
    ticks.append(len(distData))
    plt.xticks(tickloc,ticks)
    plt.yticks(tickloc,ticks)
    plt.gca().invert_yaxis()
    plt.tick_params('both', length=10, width=2)
    plt.colorbar()
    plt.show()
    os.chdir("..")
    os.chdir("..")
    os.chdir("..")



