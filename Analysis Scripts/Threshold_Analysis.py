import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import os,glob

def FindHinges(poly_size):
    xyzFile = open("../initial.xyz", "r")
    line = ""
    while line !="Angles \n":
        line = xyzFile.readline()
        if not line:
            return range(2,poly_size-1)
    line = xyzFile.readline()
    line = xyzFile.readline()
    hinges = []
    for i in range(2,poly_size):
        if line and line.split()[3] == str(i):
            line = xyzFile.readline()
            continue
        hinges.append(i)
    return hinges


def calc_rg():
    readInfo = open("Poly_Info.txt", "r")
    line = readInfo.readline()
    while line != "Time\t\tSph\t\tRg\n":
        line = readInfo.readline()
    rg = 0.0
    i = 0
    for x in readInfo:
        line_info = x.split()
        if(float(line_info[0])>1000):

            i = i + 1
            rg = rg + float(line_info[2])

    rg = rg/i
    readInfo.close()
    return rg


def CreatePAFile():
    filename = []
    for file in glob.glob("dump.*.txt"):
        filename.append(file)
    file= sorted(filename)
    total_steps=len(file)

    polymer_write = open("../AdsInfo.txt", "w+")
    a = np.array(np.genfromtxt(file[0], skip_header=9))
    polysize  = len(a)
    hinges  = FindHinges(polysize)
    for i in hinges:
        polymer_write.write(str(i)+"\t")
    polymer_write.write("\n")
    ads_Total = []

    k = 0
    for j in range(1000,total_steps,1):

        file_name=file[j]
        if(k==0):
            ads = np.zeros(polysize)
        a = np.array(np.genfromtxt(file_name, skip_header=9))
        a = np.array(sorted(np.array(a), key=lambda a_entry: a_entry[0]))        ### sort them by their monomer number

        poly_ads = ((np.abs(a[:,7])>48.5))
        ads = poly_ads+ads
        k= k+1
        if k ==100:
            ads = ads/100
            ads_Total.append(ads)
            k = 0
    #polymer_write.write("Time\t\tSph\t\tRg\n")

    ads_Total = np.array(ads_Total)
    for i in range(len(ads_Total)):
        line_str = str(ads_Total[i,0])
        for j in range(1,polysize):
            line_str = line_str+"\t"+str(ads_Total[i,j])
        polymer_write.write(" %d\t%s\n"%(i,line_str))

    polymer_write.close()

def CreatePAFilesHandle(parentFolder):
    os.chdir(os.getcwd()+"\\"+parentFolder)
    p_dirs = next(os.walk('.'))[1]
    os.chdir(os.getcwd() + "\\" + p_dirs[0])
    c_dirs = next(os.walk('.'))[1]
    os.chdir("..")

    c_dirs.sort(key = lambda x: int(x.split("_")[2]))



    for dirs_p in p_dirs:


        os.chdir(os.getcwd() + "\\" + dirs_p)
        c_dirs = next(os.walk('.'))[1]
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            print(parentFolder+"\\"+dirs_p+"\\"+dirs)
            os.chdir(os.getcwd() + "\\" + dirs+"\\coord_dump")
            CreateTotalPAFile()



            os.chdir("..")
            os.chdir("..")

        os.chdir("..")
    os.chdir("..")

def CreateTotalPAFile():
    filename = []
    for file in glob.glob("dump.*.txt"):
        filename.append(file)
    file= sorted(filename)
    total_steps=len(file)

    polymer_write = open("../AdsInfoTotal.txt", "w+")
    a = np.array(np.genfromtxt(file[0], skip_header=9))
    polysize  = len(a)
    hinges  = FindHinges(polysize)
    for i in hinges:
        polymer_write.write(str(i)+"\t")
    polymer_write.write("\n")
    ads_Total = []

    k = 0
    for j in range(1000,total_steps,1):

        file_name=file[j]
        if(k==0):
            ads = np.zeros(polysize)
        a = np.array(np.genfromtxt(file_name, skip_header=9))
        a = np.array(sorted(np.array(a), key=lambda a_entry: a_entry[0]))        ### sort them by their monomer number

        poly_ads = ((np.abs(a[:,7])>48.5))
        ads = poly_ads+ads
        k= k+1
    #polymer_write.write("Time\t\tSph\t\tRg\n")
    ads = ads/k
    ads_Total = np.array(ads_Total)
    line_str = ""
    for j in range(0,polysize):
        line_str = line_str+"\t"+str(ads[j])
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
            CreatePAFile()
            print(os.getcwd())
            os.chdir("..")
            print(os.getcwd())
            os.chdir("..")
            print(os.getcwd())
        os.chdir("..")
    os.chdir("..")

def CreatePAGraph(PolymerLocation):
    os.chdir(os.getcwd()+"\\"+PolymerLocation)
    polymer_write = open("AdsInfo.txt", "r")
    hinges = polymer_write.readline()
    hinges = hinges.split()



    data = np.array(np.genfromtxt("AdsInfo.txt", skip_header=1))
    data = data[:,1:]
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    if(len(hinges)>0):
        hinges = np.array(hinges).astype(int)
        ax.scatter3D(hinges,np.zeros(len(hinges)),np.zeros(len(hinges)),c = "red")
    #x = []
    #y = []
    #z = []
    #for i in range(len(data)):
    #    for j in range(len(data[0])):
    #        x.append(i)
    #        y.append(j)
    #        z.append(data[i,j])
    #ax.scatter3D(x,y,z,c = z, alpha = 1)

    for i in range(len(data)):
        x = range(len(data[0]))
        y = np.ones(len(data[0]))*i
        z = data[i,:]
        ax.plot3D(x,y,z)
    fig.show()
    os.chdir("..")
    os.chdir("..")
    os.chdir("..")

def CreatePATotalGraphs(parentFolder):
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

    hingeLocs = []
    polyAds = []
    init = 0
    for dirs in c_dirs:
        labels.append(int(dirs[5:]))


    for dirs_p in p_dirs:


        os.chdir(os.getcwd() + "\\" + dirs_p)
        c_dirs = next(os.walk('.'))[1]
        pAds = []
        rHinge = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)
            readInfo = open("AdsInfoTotal.txt", "r")
            hinges = readInfo.readline()
            if init ==0:

                hinges = np.array(hinges.split())
                hinges = hinges.astype(int)
                hingeLocs.append(hinges)
            ads = np.array(readInfo.readline().split())
            ads = ads.astype(np.float)

            pAds.append(ads)
            os.chdir("..")
        os.chdir("..")
        init = 1
        polyAds.append(np.array(pAds))



    cols = int(np.ceil(len(labels)/3))
    fig, axs = plt.subplots(3, cols)
    fig.suptitle("Fractional Adsorption")
    polyAds = np.array(polyAds)
    index = 0
    for i in range(0,len(polyAds[0])):
        for j in range(0,len(polyAds)):

            axs[int(np.floor((i)/cols)),int((i)%cols)].plot(range(len(polyAds[j,i])), polyAds[j,i], label = "0."+str(j) )
            axs[int(np.floor((i)/cols)),(i)%cols].legend()
            axs[int(np.floor((i)/cols)),(i)%cols].title.set_text(labels[i])
            axs[int(np.floor((i)/cols)),(i)%cols].set_ylim(0,1)
            index = index +1
        axs[int(np.floor((i)/cols)),(i)%cols].scatter(hingeLocs[i],np.zeros(len(hingeLocs[i])),s=100,c="r")
        index = 0

    fig.show()
    os.chdir("..")

def CreatePATotalGraph(parentFolder):
    os.chdir(os.getcwd()+"\\"+parentFolder)
    plt.rc('font',size = 17)
    fig, axs = plt.subplots()

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

    hingeLocs = []
    polyAds = []
    init = 0
    for dirs in c_dirs:
        labels.append(int(dirs[5:]))


    for dirs_p in p_dirs:


        os.chdir(os.getcwd() + "\\" + dirs_p)
        c_dirs = next(os.walk('.'))[1]
        pAds = []
        rHinge = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)
            print(dirs_p)
            print(dirs)
            readInfo = open("AdsInfoTotal.txt", "r")
            hinges = readInfo.readline()
            if init ==0:

                hinges = np.array(hinges.split())
                hinges = hinges.astype(int)
                hingeLocs.append(hinges)
            ads = np.array(readInfo.readline().split())
            ads = ads.astype(np.float)

            pAds.append(ads)
            os.chdir("..")
        os.chdir("..")
        init = 1
        polyAds.append(np.array(pAds))



    cols = int(np.ceil(len(labels)/3))
    #fig, axs = plt.subplots(3, cols)
    #fig.suptitle("Fractional Adsorption")
    polyAds = np.array(polyAds)
    index = 0
    value = 10
    for i in range(0,len(polyAds[0])):
        if(labels[i]==value):
            for j in range(0,len(polyAds)):

                axs.plot(range(1,1+len(polyAds[j,i])), polyAds[j,i], label = r'$\pi_c$ = 0.'+str(j), linewidth = 3)
                #axs.legend(fontsize=14, loc='center right',frameon=False, bbox_to_anchor=(1.3,.5))
                #plt.title(labels[i])
                #leg = ax2.legend(frameon=False,fontsize=14)
                axs.set_ylim(0,1)
                axs.set_xlim(1,50)
                axs.set_xlabel('Bead Index', fontsize=20)
                axs.set_ylabel(r'Fraction Adsorption ($F_A$)', fontsize=20,)
            index = index +1
            #for k in range(0,200,1):
            #for p in range(1,51):
                #plt.scatter(p,.33)
            #    axs.scatter(np.array(hingeLocs[i])-1,np.zeros(len(hingeLocs[i]))+.005*k,s=15,c="grey",alpha=.2)
            for hinge in hingeLocs[i]:
                ret = patches.Rectangle((hinge-1.5,0),2,1, facecolor="lightgrey")
                axs.add_patch(ret)
            index = 0


    plt.show()

    #fig.show()
    os.chdir("..")

def CompareStiffRegions(parentFolders):

    FAs = []
    HingeTotal = []
    x_vals = []
    labels = []
    for parentFolder in parentFolders:
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

        hingeLocs = []
        polyAds = []
        init = 0
        for dirs in c_dirs:
            labels.append(int(dirs[5:]))


        for dirs_p in p_dirs:


            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            pAds = []
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                readInfo = open("AdsInfoTotal.txt", "r")
                hinges = readInfo.readline()
                if init ==0:

                    hinges = np.array(hinges.split())
                    hinges = hinges.astype(int)
                    hingeLocs.append(hinges)
                ads = np.array(readInfo.readline().split())
                ads = ads.astype(np.float)

                pAds.append(ads)
                os.chdir("..")
            os.chdir("..")
            init = 1
            polyAds.append(np.array(pAds))
        FAs.append(polyAds)
        HingeTotal.append(hingeLocs)
        os.chdir("..")


    stiffRegion = np.zeros((len(FAs),len(np.array(FAs[0])),len(np.array((FAs[0][0])))))
    for z in range(len(FAs)):
        for i in range(len(FAs[0])):
            for j in range(len(FAs[0][0])):
                numStiff = 0
                for k in range(len(FAs[0][0][0])-2):
                    #if(k+2 not in HingeTotal[z][j]):
                    numStiff = numStiff +1
                    stiffRegion[z,i,j] = stiffRegion[z,i,j]+ FAs[z][i][j][k]
                if numStiff>0:
                    stiffRegion[z,i,j] = stiffRegion[z,i,j]/numStiff

    #for i in range(HingeTotal)
    for i in range(len(stiffRegion[0][0])):
        plt.plot(np.arange(len(stiffRegion[0]))/10,(np.array(stiffRegion)[0,:,i]-np.array(stiffRegion[1,:,i])),label=labels[i],marker = 'o')
        plt.ylabel("FA Center - FA Clustered",fontsize=20)
        plt.xlabel(r'$\pi_c$', fontsize=20)
        plt.title("Fractional Adsorption Difference of Stiff Regions")
        plt.legend()
    plt.show()
    os.chdir("..")


def GraphAdsThreshold(parentFolders,thresh=.5):
    Thresholds = []

    for parentFolder in parentFolders:
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
        init = 0
        for dirs in c_dirs:
            labels.append(int(dirs[5:])/48)
        adsorptionThresh = np.ones(len(labels))*-1
        previous = (np.ones(len(labels))*-1).astype(float)
        p = 0
        for dirs_p in p_dirs:
            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            c = 0

            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                readInfo = open("AdsInfoTotal.txt", "r")
                hinges = readInfo.readline()
                adsorption = readInfo.readline().split()
                totalAds = np.average(np.array(adsorption).astype(float))
                if totalAds > thresh and adsorptionThresh[c] == -1:
                    adsorptionThresh[c] = x_vals[p-1] + (x_vals[p]-x_vals[p-1])* ((thresh-previous[c])/(totalAds-previous[c]))
                elif(adsorptionThresh[c]==-1):
                    previous[c] = totalAds
                c = c+1
                os.chdir("..")
            p = p+1
            os.chdir("..")
        os.chdir("..")
        Thresholds.append(adsorptionThresh)
    i = 0
    k = 0
    random = []
    total_Random = []

    while i< len(Thresholds):

        if "Random" in parentFolders[i]:
            total_Random.append(Thresholds[i])
            if k == 0:
                random  = np.array(Thresholds[i])
            else:
                random = (random*k+np.array(Thresholds[i]))/(k+1)
            del Thresholds[i]
            del parentFolders[i]
            i = i-1
            k  = k+1

        i = i+1
    i=0
    total_Random = np.array(total_Random)
    std = np.sqrt(np.var(total_Random,axis=0))


    if k!=0:
        Thresholds.append(random)
        parentFolders.append("Random_hi")

    plt.rc('font',size = 17)

    for thre in Thresholds:
        #plt.plot(labels,thre,marker = 'o')
        if parentFolders[i].split("_")[0] =="Random":
            plt.errorbar(labels,thre,yerr=std,label=parentFolders[i].split("_")[0],marker = 'o')
        else:
            plt.plot(labels,thre,label=parentFolders[i].split("_")[0],marker = 'o')
        i = i+1
    #plt.legend()
    plt.xlabel(r'$\eta$',fontsize=18)
    plt.ylabel("Adsorption Threshold",fontsize=18)
    leg = plt.legend(frameon=False,fontsize=14)
    plt.show()
    print("hi")

def GraphRgThreshold(parentFolders,thresh=.5,graphTypes = []):
    Thresholds = []
    plt.rc('font',size = 17)
    for parentFolder in parentFolders:
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
        init = 0
        Rg_0 = []
        for dirs in c_dirs:
            labels.append(int(dirs[5:])/48)

            os.chdir(os.getcwd()+"\\p_0\\"+dirs)
            Rg_0.append(calc_rg())
            os.chdir("..")
            os.chdir("..")

        RgFinal = np.ones(len(labels))*-1
        previous = (np.ones(len(labels))*-1).astype(float)
        p = 0
        """
        for dirs_p in p_dirs:
            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            c = 0

            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                Rg = calc_rg()/Rg_0[c]
                if Rg > thresh and RgFinal[c] == -1:
                    RgFinal[c] = x_vals[p-1] + (x_vals[p]-x_vals[p-1])* ((thresh-previous[c])/(Rg-previous[c]))
                elif(RgFinal[c]==-1):
                    previous[c] = Rg
                c = c+1
                os.chdir("..")
            p = p+1
            os.chdir("..")
        os.chdir("..")"""
        os.chdir(os.getcwd()+"\\p_0.9")
        c_dirs = next(os.walk('.'))[1]
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        c = 0
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)

            RgFinal[c] = calc_rg()/Rg_0[c]

            c = c+1
            os.chdir("..")

        Thresholds.append(RgFinal)
        os.chdir("..")
        os.chdir("..")

        i = 0
    k = 0
    random = []
    random_NW = []
    b = 0
    random_Total = []
    while i< len(Thresholds):

        if "Random" in parentFolders[i]:
            if "NW" in parentFolders[i]:
                random_Total.append(Thresholds[i])
                if b == 0:
                    random_NW  = np.array(Thresholds[i])
                else:
                    random_NW = (random_NW*b+np.array(Thresholds[i]))/(b+1)
                del Thresholds[i]
                del parentFolders[i]

                i = i-1
                b  = b+1
            else:
                if k == 0:
                    random  = np.array(Thresholds[i])
                else:
                    random = (random*k+np.array(Thresholds[i]))/(k+1)
                del Thresholds[i]
                del parentFolders[i]
                i = i-1
                k  = k+1
        i = i+1

    random_Total = np.array(random_Total)
    var = np.sqrt(np.var(random_Total,axis=0))
    print(np.sqrt(var))
    print(random_Total)

    if k!=0:
        Thresholds.append(random)
        parentFolders.append("Random_hi")
    if b!= 0:
        Thresholds.append(random_NW)
        parentFolders.append("Random_NW")

    i = 0

    j =0
    k = 0

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    #plt.plot(labels, Thresholds[4],marker = 'o', label= "Bulk", c = 'black' )
    #plt.plot(labels, Thresholds[4],marker = 'o',linestyle="--", label= "Surface" ,c = 'black')

    for thre in Thresholds[:-1]:
        #plt.plot(labels,thre,marker = 'o')
        if "NW" not in parentFolders[i]:
            plt.errorbar(labels,thre,yerr=var,fmt='o',marker = 'o',linestyle="--",c =  colors[j])
            j +=1
        else:

            plt.plot(labels,thre, label=parentFolders[i].split("_")[0],marker = 'o',c = colors[k])
            k +=1
        i = i+1
    plt.errorbar(labels,Thresholds[-1],yerr = var,label=parentFolders[-1].split("_")[0],marker = 'o',c = colors[k] )
    #plt.legend()
    plt.xlabel(r'$\eta$',fontsize=18)
    plt.ylabel(r'$\langle R_g(0.9)\rangle /\langle R_g(0)\rangle$', fontsize=18)

    leg = plt.legend(frameon=False,fontsize=14)
    plt.show()
    




