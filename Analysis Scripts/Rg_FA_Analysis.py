#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import os,glob
import matplotlib as mpl


def no_within(b,z_min,z_max,sigma):
    return(len(list(x for x in b if ((x>(z_max-1))|(x<(z_min+sigma))))))


def AnalizeInfo():
    filename, z, com, z_coord, monomer, rg = [], [], [], [], [], []
    for file in glob.glob("dump.*.txt"):
        filename.append(file)
    file= sorted(filename)
    total_steps=len(file)

    t=range(0,50*total_steps,50)  ##with timestep of 0.005 and "record" timesteps every 10000

    polymer_write = open("../Poly_Info.txt", "w+")
    sph = []
    a = np.array(np.genfromtxt(file[0], skip_header=9))
    polysize  = len(a)
    total_cos = np.zeros(polysize - 2)
    total_angle = np.zeros(polysize - 2)
    total_distance = np.zeros(polysize - 1)
    percent_absorbed = 0

    for j in range(0,total_steps,1):
        #print (file[j])
        file_name=file[j]

        if j>1000:
            a = np.array(np.genfromtxt(file_name, skip_header=9))
            a = np.array(sorted(np.array(a), key=lambda a_entry: a_entry[0]))        ### sort them by their monomer number
            number_wa=no_within(a[:,5],-51,51,1.0)  ##Find out the number of monomers within 1 sigma of -50 and 50
            Sij = np.zeros((3, 3))

            com.append(np.array([np.mean(a[:, m]) for m in range(5,8)]))  ## Find out the center of mass
            #test = a[:,5:8]-com[j-1001]
            rg.append(np.sum(np.sum(np.square(a[:,5:8]-com[j-1001]),axis =1))/len(a[:,0])) ## Calculate Radius of Gyration
            for n in range(len(a)):
                for i in range(3):
                    for k in range(3):
                        Sij[i,k] = Sij[i,k]+(a[n,(5+i)]-com[j-1001][i])*(a[n,5+k]-com[j-1001][k])
            Sij = Sij/len(a)
            [eig, eigvec] = np.linalg.eig(Sij)
            eig = -np.sort(-eig)
            b = eig[0] - (1/2)*(eig[1]+eig[2])
            sph.append(b)
            z.append(number_wa/51)    ## Fraction of monomers within the wall


            dist_vector = a[1:,5:8] - a[:-1,5:8]
            total_distance = total_distance+np.sqrt(np.sum(np.square(dist_vector),axis = 1))


            dot_vector = []
            mag_vector = []
            for i in dist_vector[1:]:
                dot_vector.append(np.dot(dist_vector[0],i))
                mag_vector.append(np.sqrt(np.sum(np.square(i)))*np.sqrt(np.sum(np.square(dist_vector[0]))))
        #mag_vector = np.sum(np.square(dist_vector),axis=1)
            cos_ang = np.array(dot_vector)/np.array(mag_vector)
            angle = np.arccos(cos_ang)
            total_angle = angle + total_angle
            total_cos = total_cos +cos_ang


            #holds = (np.abs(a[:,7])>48.5)
            percent_absorbed = percent_absorbed + np.sum((np.abs(a[:,7])>48.5))/polysize

    total_distance = total_distance/(total_steps-1000)
    total_cos = total_cos/(total_steps-1000)
    total_angle = total_angle/(total_steps-1000)
    percent_absorbed = percent_absorbed/(total_steps-1000)
    index = len(total_cos)
    if(np.argmax(total_cos<0)>1):
        index = np.argmax(total_cos<0)
    total_cos_no_neg = total_cos[:index]
    total_dist_no_neg = total_distance[:index]
    x = np.zeros(len(total_dist_no_neg))
    for i in range(1,len(x)):
        x[i] = x[i-1]+total_dist_no_neg[i-1]
    fit = np.polyfit(x,np.log(total_cos_no_neg),1)
    lp = -np.average(total_dist_no_neg)/fit[0]

    polymer_write.write("pl\t%.3f\n\n"%(lp))
    polymer_write.write("abs\t%.3f\n\n"%(percent_absorbed))
    polymer_write.write("Time\t\tSph\t\tRg\n")

    for i in range(total_steps-1001):
        polymer_write.write(" %d\t%8.3f\t%8.3f\n"%(i,sph[i],rg[i]))
    polymer_write.close()

def createInfo(parentFolder):
    print(os.getcwd())
    os.chdir(os.getcwd()+"\\"+parentFolder)
    print(os.getcwd())
    p_dirs = next(os.walk('.'))[1]
    print(p_dirs)
    for folder in p_dirs:
        os.chdir(os.getcwd()+"\\"+folder)
        i_dirs = next(os.walk('.'))[1]
        for i_folder in i_dirs:
            os.chdir(os.getcwd()+"\\"+i_folder+"\\coord_dump")
            AnalizeInfo()
            print(os.getcwd())
            os.chdir("..")
            print(os.getcwd())
            os.chdir("..")
            print(os.getcwd())
        os.chdir("..")
    os.chdir("..")

def CreateRgGraphs(parentFolder):
    os.chdir(os.getcwd()+"\\"+parentFolder)
    plt.rc('font',size = 17)

    p_dirs = next(os.walk('.'))[1]
    x_vals = []
    for dirs in p_dirs:
        x_vals.append(float(dirs[2:]))
    os.chdir(os.getcwd() + "\\" + p_dirs[0])
    c_dirs = next(os.walk('.'))[1]
    labels = []
    c_dirs.sort(key = lambda x: int(x.split("_")[2]))
    for dirs in c_dirs:
        labels.append(int(dirs[5:]))
    sph_0 = []
    rg_0 = []
    sph_t = []
    rg_t = []
    c_dirs.sort(key = lambda x: int(x.split("_")[2]))
    for dirs in c_dirs:
        os.chdir(os.getcwd() + "\\" + dirs)

        [sph_min, sph_max, sph, rg_min, rg_max, rg] = calc_sph_rg()
        sph_0.append(sph)
        rg_0.append(rg)
        os.chdir("..")
    os.chdir("..")
    sph_0 = np.array(sph_0)
    rg_0 = np.array((rg_0))

    for dirs_p in p_dirs:
        sph_temp = []
        rg_temp = []

        os.chdir(os.getcwd() + "\\" + dirs_p)
        c_dirs = next(os.walk('.'))[1]
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)
            [sph_min, sph_max, sph, rg_min, rg_max, rg] = calc_sph_rg()

            sph_temp.append(sph)
            rg_temp.append(rg)
            os.chdir("..")
        os.chdir("..")
        sph_t.append(np.array(sph_temp)/sph_0)
        rg_t.append(np.array(rg_temp)/rg_0)
    fig, (ax2) = plt.subplots(1)
    sph_t = np.array(sph_t)
    rg_t = np.array(rg_t)

   # for i in range(len(sph_t[0])):
   #     ax1.plot(x_vals,sph_t[:,i],label = labels[i])

    #colors = plt.cm.winter(np.linspace(0,1,len(rg_t[0])))
    colors = []
    for i in range(len(rg_t[0])):
         p = ax2.plot(x_vals,rg_t[:,i],label = r'$\eta$ = '+str(round(float(labels[i])/48,2)), marker = 'o')
         colors.append(p[0].get_color())
    #ax2.legend(fontsize=12,loc = 'lower left', title="Hinge Count")
    ax2.set_xlabel(r'Osmotic Pressure, $\Pi_c$', fontsize=18)
    ax2.set_ylabel(r'$\langle R_g(\Pi_c)\rangle /\langle R_g(0)\rangle$', fontsize=18)
    #leg = ax2.legend(frameon=False,fontsize=14)
    #for color,text in zip(colors, leg.get_texts()):
        #text.set_color(color)
    fig.show()
    os.chdir("..")

def CreatePAGraphs(parentFolder):
    os.chdir(os.getcwd()+"\\"+parentFolder)
    plt.rc('font',size = 17)

    p_dirs = next(os.walk('.'))[1]
    x_vals = []
    for dirs in p_dirs:
        x_vals.append(float(dirs[2:]))
    x_vals.sort()
    os.chdir(os.getcwd() + "\\" + p_dirs[0])
    c_dirs = next(os.walk('.'))[1]
    os.chdir("..")

    labels = []
    c_dirs.sort(key = lambda x: int(x.split("_")[2]))

    for dirs in c_dirs:
        labels.append(int(dirs[5:]))


    pAbs_t = []


    for dirs_p in p_dirs:


        os.chdir(os.getcwd() + "\\" + dirs_p)
        c_dirs = next(os.walk('.'))[1]
        pAbs = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)
            readInfo = open("Poly_Info.txt", "r")
            line = readInfo.readline()
            while "abs" not in line:
                line = readInfo.readline()

            line_info = line.split()
            pAbs.append(float(line_info[1]))
            os.chdir("..")
        os.chdir("..")
        pAbs_t.append(np.array(pAbs))

    fig, ax1= plt.subplots(1)

    pAbs_t = np.array(pAbs_t)

    #test = abst[0,0]
    colors = plt.cm.winter(np.linspace(0,1,len(pAbs_t[0])))
    for i in range(len(pAbs_t[0])):
        ax1.plot(x_vals, pAbs_t[:,i],label = labels[i], marker = 'o')

    #ax1.legend(fontsize=15, title = "Hinge Count")

    ax1.set_xlabel(r'Osmotic Pressure, $\Pi_c$', fontsize=18)
    ax1.set_ylabel(r'$\langle N_{wall}/N \rangle$', fontsize=20)
    fig.show()
    os.chdir("..")

def CreateAverageRgGraphs(parentFolders):
    count = 0
    rg_compiled = []
    for parentFolder in parentFolders:
        count = count +1
        os.chdir(os.getcwd()+"\\"+parentFolder)
        p_dirs = next(os.walk('.'))[1]
        x_vals = []
        for dirs in p_dirs:
            x_vals.append(float(dirs[2:]))
        os.chdir(os.getcwd() + "\\" + p_dirs[0])
        c_dirs = next(os.walk('.'))[1]
        labels = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            labels.append(int(dirs[5:]))
        sph_0 = []
        rg_0 = []
        sph_t = []
        rg_t = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)

            [sph_min, sph_max, sph, rg_min, rg_max, rg] = calc_sph_rg()
            sph_0.append(sph)
            rg_0.append(rg)
            os.chdir("..")
        os.chdir("..")
        sph_0 = np.array(sph_0)
        rg_0 = np.array((rg_0))

        for dirs_p in p_dirs:
            sph_temp = []
            rg_temp = []

            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                [sph_min, sph_max, sph, rg_min, rg_max, rg] = calc_sph_rg()

                sph_temp.append(sph)
                rg_temp.append(rg)
                os.chdir("..")
            os.chdir("..")
            sph_t.append(np.array(sph_temp)/sph_0)
            rg_t.append(np.array(rg_temp)/rg_0)
        rg_compiled.append(rg_t)
        os.chdir("..")

    RgMin = np.zeros((len(rg_compiled[0]),len(rg_compiled[0][0])))+50
    RgMax = np.zeros((len(rg_compiled[0]),len(rg_compiled[0][0])))

    pRg = np.array(rg_compiled[0])/count
    for i in range(1,len(rg_compiled)):
        for j in range(len(rg_compiled[0])):
            for k in range(len(rg_compiled[0][0])):
                if i !=0:
                    pRg[j,k]= rg_compiled[i][j][k]/count+pRg[j,k]
                if rg_compiled[i][j][k]> RgMax[j,k]:
                    RgMax[j][k] = rg_compiled[i][j][k]
                if rg_compiled[i][j][k] < RgMin[j,k]:
                    RgMin[j][k] = rg_compiled[i][j][k]


    fig, (ax2) = plt.subplots(1)
    sph_t = np.array(sph_t)
    rg_t = np.array(rg_t)
    # for i in range(len(sph_t[0])):
    #     ax1.plot(x_vals,sph_t[:,i],label = labels[i])
    avgDiff = (RgMax-RgMin)/2
    for i in range(len(rg_t[0])):
        if i != -5 and i!=-1:
            ax2.errorbar(x_vals,pRg[:,i],label=labels[i], marker='o',capsize=2,color='C'+str(i))
            #ax2.scatter(x_vals,RgMin[:,i],marker = '*', color='C'+str(i))
            #ax2.scatter(x_vals,RgMax[:,i],marker = '*', color= 'C'+str(i))
    #ax2.legend(fontsize=12,loc = 'lower left', title="Hinge Count")
    ax2.set_xlabel(r'$\pi_c$', fontsize=18)
    ax2.set_ylabel("Normalized Radius of Gyration", fontsize=18)
    fig.show()


def CreateRgHist(parentFolder):

    rg_compiled = []
    os.chdir(os.getcwd()+"\\"+parentFolder)
    p_dirs = next(os.walk('.'))[1]
    x_vals = []
    for dirs in p_dirs:
        x_vals.append(float(dirs[2:]))
    os.chdir(os.getcwd() + "\\" + p_dirs[0])
    c_dirs = next(os.walk('.'))[1]
    labels = []
    c_dirs.sort(key = lambda x: int(x.split("_")[2]))
    for dirs in c_dirs:
        labels.append(int(dirs[5:]))
    average_rg_0 = []

    freq_t = []
    c_dirs.sort(key = lambda x: int(x.split("_")[2]))
    for dirs in c_dirs:
        os.chdir(os.getcwd() + "\\" + dirs)

        [_, _, _, _, _, rg] = calc_sph_rg()
        average_rg_0.append(rg)
        os.chdir("..")
    os.chdir("..")

    average_rg_0 = np.array((average_rg_0))

    for dirs_p in p_dirs:
        freq_temp = []
        os.chdir(os.getcwd() + "\\" + dirs_p)
        c_dirs = next(os.walk('.'))[1]
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))

        for dirs in c_dirs:
            freq  = np.zeros((300))
            os.chdir(os.getcwd() + "\\" + dirs)
            rg = np.array(GetRg())/average_rg_0[len(freq_temp)]
            for value in rg:
                freq[int(round(value*100)/5)] +=1
            freq_temp.append(freq)
            os.chdir("..")
        os.chdir("..")
        #sph_t.append(np.array(sph_temp)/sph_0)
        freq_t.append(np.array(freq_temp)/average_rg_0[:,None])




    os.chdir("..")

    freq_t = np.array(freq_t)
    cols = int(np.ceil(10/3))
    fig, axs = plt.subplots(3, cols)
    for i, p_value in enumerate(freq_t):
        for j, freq_v in enumerate(p_value):


            axs[int(np.floor((i)/cols)),int((i)%cols)].plot(5*np.arange(len(freq_v))/100,freq_v/sum(freq_v), label = labels[j] )
            axs[int(np.floor((i)/cols)),(i)%cols].legend()
            axs[int(np.floor((i)/cols)),(i)%cols].title.set_text("0."+str(i))
            axs[int(np.floor((i)/cols)),(i)%cols].set_xlim(0,2)

        #axs[int(np.floor((i)/cols)),(i)%cols].scatter(hingeLocs[i],np.zeros(len(hingeLocs[i])),s=100,c="r")




    """
    pRg = np.array(rg_compiled[0])/count
    for i in range(1,len(rg_compiled)):
        for j in range(len(rg_compiled[0])):
            for k in range(len(rg_compiled[0][0])):
                if i !=0:
                    pRg[j,k]= rg_compiled[i][j][k]/count+pRg[j,k]
                if rg_compiled[i][j][k]> RgMax[j,k]:
                    RgMax[j][k] = rg_compiled[i][j][k]
                if rg_compiled[i][j][k] < RgMin[j,k]:
                    RgMin[j][k] = rg_compiled[i][j][k]


    fig, (ax2) = plt.subplots(1)
    sph_t = np.array(sph_t)
    freq_t = np.array(freq_t)
    # for i in range(len(sph_t[0])):
    #     ax1.plot(x_vals,sph_t[:,i],label = labels[i])
    avgDiff = (RgMax-RgMin)/2
    for i in range(len(freq_t[0])):
        if i != -5 and i!=-1:
            ax2.errorbar(x_vals,pRg[:,i],label=labels[i], marker='o',capsize=2,color='C'+str(i))
            #ax2.scatter(x_vals,RgMin[:,i],marker = '*', color='C'+str(i))
            #ax2.scatter(x_vals,RgMax[:,i],marker = '*', color= 'C'+str(i))
    #ax2.legend(fontsize=12,loc = 'lower left', title="Hinge Count")
    ax2.set_xlabel(r'$\pi_c$', fontsize=18)
    ax2.set_ylabel("Normalized Radius of Gyration", fontsize=18)
    fig.show()"""




def CreateAveragePAGraphs(parentFolders):
    count = 0
    abs_Compiled = []
    for parentFolder in parentFolders:
        count = count + 1
        os.chdir(os.getcwd()+"\\"+parentFolder)
        p_dirs = next(os.walk('.'))[1]
        x_vals = []
        for dirs in p_dirs:
            x_vals.append(float(dirs[2:]))
        x_vals.sort()
        os.chdir(os.getcwd() + "\\" + p_dirs[0])
        c_dirs = next(os.walk('.'))[1]
        os.chdir("..")

        labels = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))

        for dirs in c_dirs:
            labels.append(int(dirs[5:]))


        pAbs_t = []


        for dirs_p in p_dirs:
            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            pAbs = []
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                readInfo = open("Poly_Info.txt", "r")
                line = readInfo.readline()
                while "abs" not in line:
                    line = readInfo.readline()

                line_info = line.split()
                pAbs.append(float(line_info[1]))
                os.chdir("..")
            os.chdir("..")
            pAbs_t.append(np.array(pAbs))
        abs_Compiled.append(pAbs_t)
        os.chdir("..")
    fig, ax1= plt.subplots(1)


    pAbs = np.array(abs_Compiled[0])/count
    AbsMin = np.zeros((len(abs_Compiled[0]),len(abs_Compiled[0][0])))+50
    AbsMax = np.zeros((len(abs_Compiled[0]),len(abs_Compiled[0][0])))


    for i in range(len(abs_Compiled)):
        for j in range(len(abs_Compiled[0])):
            for k in range(len(abs_Compiled[0][0])):
                if i !=0:
                    pAbs[j,k]= abs_Compiled[i][j][k]/count+pAbs[j,k]
                if abs_Compiled[i][j][k]> AbsMax[j,k]:
                    AbsMax[j][k] = abs_Compiled[i][j][k]
                if abs_Compiled[i][j][k] < AbsMin[j,k]:
                    AbsMin[j][k] = abs_Compiled[i][j][k]

    avgDiff = (AbsMax-AbsMin)/2
    #test = abst[0,0]
    #for i in range(len(pAbs[0])):
        #ax1.plot(x_vals, pAbs[:,i],label = labels[i], marker = 'o')
    #    if i != -5 and i!=-1:
    #        ax1.errorbar(x_vals,pAbs[:,i],label=labels[i], marker='o',capsize=2,color='C'+str(i))
            #ax1.scatter(x_vals,AbsMin[:,i],marker = '*', color='C'+str(i))
            #ax1.scatter(x_vals,AbsMax[:,i],marker = '*', color= 'C'+str(i))

    #ax1.legend(fontsize=15, title = "Hinge Count")

    ax1.set_xlabel(r'$\pi_c$', fontsize=20)
    ax1.set_ylabel("Fraction Adsorption", fontsize=20)
    fig.show()

def CreateSubplotRgGraphs(parentFolders,fig,axs):
    count = 0
    rg_compiled = []
    for parentFolder in parentFolders:
        count = count +1
        os.chdir(os.getcwd()+"\\"+parentFolder)
        p_dirs = next(os.walk('.'))[1]
        x_vals = []
        for dirs in p_dirs:
            x_vals.append(float(dirs[2:]))
        os.chdir(os.getcwd() + "\\" + p_dirs[0])
        c_dirs = next(os.walk('.'))[1]
        labels = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            labels.append(int(dirs[5:]))
        sph_0 = []
        rg_0 = []
        sph_t = []
        rg_t = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)

            [sph_min, sph_max, sph, rg_min, rg_max, rg] = calc_sph_rg()
            sph_0.append(sph)
            rg_0.append(rg)
            os.chdir("..")
        os.chdir("..")
        sph_0 = np.array(sph_0)
        rg_0 = np.array((rg_0))

        for dirs_p in p_dirs:
            sph_temp = []
            rg_temp = []

            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                [sph_min, sph_max, sph, rg_min, rg_max, rg] = calc_sph_rg()

                sph_temp.append(sph)
                rg_temp.append(rg)
                os.chdir("..")
            os.chdir("..")
            sph_t.append(np.array(sph_temp)/sph_0)
            rg_t.append(np.array(rg_temp)/rg_0)
        rg_compiled.append(rg_t)
        os.chdir("..")

    RgMin = np.zeros((len(rg_compiled[0]),len(rg_compiled[0][0])))+50
    RgMax = np.zeros((len(rg_compiled[0]),len(rg_compiled[0][0])))

    pRg = np.array(rg_compiled[0])/count
    for i in range(1,len(rg_compiled)):
        for j in range(len(rg_compiled[0])):
            for k in range(len(rg_compiled[0][0])):
                if i !=0:
                    pRg[j,k]= rg_compiled[i][j][k]/count+pRg[j,k]
                if rg_compiled[i][j][k]> RgMax[j,k]:
                    RgMax[j][k] = rg_compiled[i][j][k]
                if rg_compiled[i][j][k] < RgMin[j,k]:
                    RgMin[j][k] = rg_compiled[i][j][k]



    cols = int(np.ceil(len(labels)/3))
    #fig, axs = plt.subplots(1, 5)
    #fig.suptitle("Radius of Gyration")
    index = 0
    rg_compiled = np.array(rg_compiled)
    for i in range(1,len(rg_compiled[0,0])-1):
        for j in range(0,len(rg_compiled)):
            if i!=1:
                axs[0,int(i)-1].yaxis.set_tick_params(labelleft=False)
                axs[0,int(i)-1].set_yticks([])
            else:
                axs[0,int(i)-1].yaxis.set_tick_params(labelsize = 17)

            axs[0,int(i)-1].xaxis.set_tick_params(labelleft=False)
            axs[0,int(i)-1].set_xticks([])

            axs[0,int(i)-1].plot(x_vals, rg_compiled[j,:,i], label = "Run "+str(j) )
            #axs[int(np.floor((i)/cols)),(i)%cols].legend()
            #axs[int(np.floor((i)/cols)),(i)%cols].title.set_text(labels[i])
            axs[0,int(i)-1].set_ylim(0,1.2)
            index = index +1

        index = 0
    #fig.show()

def CreateSubplotPAGraphs(parentFolders,fig,axs):
    count = 0
    abs_Compiled = []
    for parentFolder in parentFolders:
        count = count + 1
        os.chdir(os.getcwd()+"\\"+parentFolder)
        p_dirs = next(os.walk('.'))[1]
        x_vals = []
        for dirs in p_dirs:
            x_vals.append(float(dirs[2:]))
        x_vals.sort()
        os.chdir(os.getcwd() + "\\" + p_dirs[0])
        c_dirs = next(os.walk('.'))[1]
        os.chdir("..")

        labels = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))

        for dirs in c_dirs:
            labels.append(int(dirs[5:]))


        pAbs_t = []


        for dirs_p in p_dirs:
            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            pAbs = []
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                readInfo = open("Poly_Info.txt", "r")
                line = readInfo.readline()
                while "abs" not in line:
                    line = readInfo.readline()

                line_info = line.split()
                pAbs.append(float(line_info[1]))
                os.chdir("..")
            os.chdir("..")
            pAbs_t.append(np.array(pAbs))
        abs_Compiled.append(pAbs_t)
        os.chdir("..")
    #fig, ax1= plt.subplots(1)


    pAbs = np.array(abs_Compiled[0])/count
    AbsMin = np.zeros((len(abs_Compiled[0]),len(abs_Compiled[0][0])))+50
    AbsMax = np.zeros((len(abs_Compiled[0]),len(abs_Compiled[0][0])))


    for i in range(len(abs_Compiled)):
        for j in range(len(abs_Compiled[0])):
            for k in range(len(abs_Compiled[0][0])):
                if i !=0:
                    pAbs[j,k]= abs_Compiled[i][j][k]/count+pAbs[j,k]
                if abs_Compiled[i][j][k]> AbsMax[j,k]:
                    AbsMax[j][k] = abs_Compiled[i][j][k]
                if abs_Compiled[i][j][k] < AbsMin[j,k]:
                    AbsMin[j][k] = abs_Compiled[i][j][k]

    cols = int(np.ceil(len(labels)/3))
    #fig, axs = plt.subplots(1,5)
    fig.suptitle("Fraction Adsorption")
    index = 0
    abs_Compiled = np.array(abs_Compiled)
    plt.subplots_adjust(wspace=.1,hspace = .1)

    for i in range(1,len(abs_Compiled[0,0])-1):
        for j in range(0,len(abs_Compiled)):
            if i!=1:
                axs[1,int(i)-1].yaxis.set_tick_params(labelleft=False)
                axs[1,int(i)-1].set_yticks([])
            axs[1,int(i)-1].yaxis.set_tick_params(labelsize = 17)
            axs[1,int(i)-1].xaxis.set_tick_params(labelsize = 17)

            axs[1,int(i)-1].plot(x_vals, abs_Compiled[j,:,i], label = "Run "+str(j) )
            axs[1,int(i)-1].set_xticks([0,0.2,0.4,0.6,0.8])
            #axs[int(np.floor((i)/cols)),(i)%cols].legend()
            #axs[int(np.floor((i)/cols)),(i)%cols].title.set_text(labels[i])
            axs[1,int(i)-1].set_ylim(0,1.2)
            index = index +1

        index = 0
    #fig.show()

def calc_sph_rg():
    readInfo = open("Poly_Info.txt", "r")
    line = readInfo.readline()
    while line != "Time\t\tSph\t\tRg\n":
        line = readInfo.readline()
    sph = 0.0
    sph_min = 999
    sph_max = 0

    rg_min = 999
    rg_max = 0
    rg = 0.0
    i = 0
    for x in readInfo:
        line_info = x.split()
        if(float(line_info[0])>1):

            sph_min = min(float(line_info[1]),sph_min)
            sph_max = max(float(line_info[1]),sph_max)

            rg_min = min(float(line_info[2]),rg_min)
            rg_max = max(float(line_info[2]),rg_max)

            i = i + 1
            sph = sph + float(line_info[1])
            rg = rg + float(line_info[2])
    sph = sph/i
    rg = rg/i
    readInfo.close()
    return [sph_min, sph_max, sph, rg_min, rg_max, rg]

def GetRg():
    readInfo = open("Poly_Info.txt", "r")
    line = readInfo.readline()
    while line != "Time\t\tSph\t\tRg\n":
        line = readInfo.readline()



    rg  = []
    i = 0
    for x in readInfo:
        line_info = x.split()
        if(float(line_info[0])>1):
            i = i + 1
            rg.append(float(line_info[2]))
    readInfo.close()
    return rg

def calc_rand_var(parentFolder):
    os.chdir(os.getcwd() + "\\" + parentFolder)
    p_dirs = next(os.walk('.'))[1]
    rand_var_t = []
    num_atoms = 100
    hinge_num_list = []
    pic_list = p_dirs

    for dirs_p in p_dirs:
        os.chdir(os.getcwd() + "\\" + dirs_p)
        c_dirs = next(os.walk('.'))[1]
        rand_var = []
        for dir_c in c_dirs:
            hinge_num_list = c_dirs

            os.chdir(os.getcwd() + "\\" + dir_c)
            info = open("initial.xyz")
            line = info.readline()
            hingeLoc = []
            while "Angles" not in line:
                if(line==''):
                    rand_var.append(0)
                    info.close()
                    os.chdir("..")
                    break
                if("atoms" in line):
                    num_atoms = int(line.split()[0])
                line = info.readline()
            if(line == ''):
                break
            line = info.readline()
            found = False
            for i in range(1, num_atoms-1):
                if( not found):
                    line = info.readline()
                found = False
                if "1" not in line:
                    break
                if(line.split()[2] != str(i)):
                    hingeLoc.append(i)
                    found = True

            if(len(hingeLoc) == 0):
                rand_var.append(0)
                os.chdir("..")
                continue
            hingeLoc = np.array(hingeLoc)

            rand_var.append(hingeLoc.var())
            info.close()
            os.chdir("..")
            print(os.getcwd())
        rand_var_t.append(rand_var)
        os.chdir("..")
    rand_var_t = np.round(rand_var_t)
    print("\t"+makeString(hinge_num_list))
    for i in range(len(rand_var_t)):
        print(pic_list[i]+" "+makeString(rand_var_t[i]))
    os.chdir("..")

def makeString(list):
    stringf = ""
    for i in (list):
        stringf = stringf+" " +str(i)
    return stringf



