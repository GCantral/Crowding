
import numpy as np
import matplotlib.pyplot as plt
import os,glob
import math
import matplotlib as mpl


def get_rg():
    readInfo = open("Poly_Info.txt", "r")
    line = readInfo.readline()
    while line != "Time\t\tSph\t\tRg\n":
        line = readInfo.readline()
    rg_max = 0
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

def normalizeHeatmaps(heatmap,percentages):
    newHM = np.zeros((len(heatmap),100))
    for i in range(0,len(heatmap)):
        for j in range(0,100):
            newHM[i,j] = heatmap[i, np.argmin(np.abs(percentages-j/100))]
    return newHM

def createRgGraphs(parentFolders):
    rg_all = []
    labels = []
    for parentFolder in parentFolders:
        os.chdir(os.getcwd()+"\\"+parentFolder)
        p_dirs = next(os.walk('.'))[1]
        x_vals = []
        for dirs in p_dirs:
            x_vals.append(float(dirs[2:]))
        os.chdir(os.getcwd() + "\\" + p_dirs[0])
        c_dirs = next(os.walk('.'))[1]
        label = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            label.append(int(dirs[5:]))
        labels.append(label)
        rg_0 = []
        rg_t = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)

            rg = get_rg()
            rg_0.append(rg)
            os.chdir("..")
        os.chdir("..")
        rg_0 = np.array((rg_0))

        for dirs_p in p_dirs:
            rg_temp = []

            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                rg = get_rg()
                rg_temp.append(rg)
                os.chdir("..")
            os.chdir("..")
            rg_t.append(np.array(rg_temp)/rg_0)
        rg_t = np.array(rg_t)
        rg_all.append(rg_t)
        os.chdir("..")

    print(labels)

    fig, axs = plt.subplots(3, 3)
    index = 0
    sizes = []

    #fig.suptitle("Radius of Gyration")


    for dirs in parentFolders:
        sizes.append((dirs.split("_")[0]))

    for i in range(len(labels)):
        labels[i] = np.array(labels[i])/labels[i][len(labels[i])-1]
    index = 0
    Random = []
    rg_plot = []
    label_g = []
    for i in range(len(labels)):
        if "Random" in sizes[i]:
            Random.append(rg_all[i])
        else:
            rg_plot.append(rg_all[i])
            label_g.append(sizes[i])
    if Random:
        rg_plot.append(np.mean(np.array(Random),axis=0))
        label_g.append("Random")


    for i in range(1,len(rg_all[0])):
        for array in rg_plot:
            axs[math.floor((i-1)/3),(i-1)%3].plot(labels[index], array[i], label = label_g[index], marker="o" )
            #axs[math.floor((i-1)/3),(i-1)%3].legend()
            axs[math.floor((i-1)/3),(i-1)%3].title.set_text("0."+str(i))
            axs[math.floor((i-1)/3),(i-1)%3].set_ylim(0,1.2)
            index = index +1
        index = 0

    fig.show()


def createPAGraphs(parentFolders):
    pa_all = []
    labels = []
    for parentFolder in parentFolders:
        os.chdir(os.getcwd()+"\\"+parentFolder)
        p_dirs = next(os.walk('.'))[1]
        x_vals = []
        for dirs in p_dirs:
            x_vals.append(float(dirs[2:]))
        os.chdir(os.getcwd() + "\\" + p_dirs[0])
        c_dirs = next(os.walk('.'))[1]
        label = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            label.append(int(dirs[5:]))
        labels.append(label)
        rg_0 = []
        pa_t = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)

            rg = get_rg()
            rg_0.append(rg)
            os.chdir("..")
        os.chdir("..")
        for dirs_p in p_dirs:
            rg_temp = []

            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)

                readInfo = open("Poly_Info.txt", "r")
                line = readInfo.readline()
                while "abs" not in line:
                    line = readInfo.readline()

                line_info = line.split()
                rg_temp.append(float(line_info[1]))
                os.chdir("..")
            os.chdir("..")
            pa_t.append(np.array(rg_temp))
        pa_t = np.array(pa_t)
        pa_all.append(pa_t)
        os.chdir("..")

    print(labels)

    fig, axs = plt.subplots(3, 3)
    index = 0
    sizes = []
    #fig.suptitle("Fractional Adsorption")


    for dirs in parentFolders:
        sizes.append((dirs.split("_")[0]))

        Random = []
    pa_plot = []
    label_g = []
    for i in range(len(labels)):
        if "Random" in sizes[i]:
            Random.append(pa_all[i])
        else:
            pa_plot.append(pa_all[i])
            label_g.append(sizes[i])
    if Random:
        pa_plot.append(np.mean(np.array(Random),axis=0))
        label_g.append("Random")


    for i in range(len(labels)):
        labels[i] = np.array(labels[i])/labels[i][len(labels[i])-1]
    index = 0
    for i in range(1,len(pa_all[0])):
        for array in pa_plot:
            axs[math.floor((i-1)/3),(i-1)%3].plot(labels[index], array[i], label = label_g[index] , marker="o")
            #axs[math.floor((i-1)/3),(i-1)%3].legend()
            axs[math.floor((i-1)/3),(i-1)%3].title.set_text("0."+str(i))
            axs[math.floor((i-1)/3),(i-1)%3].set_ylim(0,1)
            index = index +1
        index = 0

    fig.show()


def createPAHeatmap(parentFolders):
    pa_all = []
    labels = []
    for parentFolder in parentFolders:
        os.chdir(os.getcwd()+"\\"+parentFolder)
        p_dirs = next(os.walk('.'))[1]
        x_vals = []
        for dirs in p_dirs:
            x_vals.append(float(dirs[2:]))
        os.chdir(os.getcwd() + "\\" + p_dirs[0])
        c_dirs = next(os.walk('.'))[1]
        label = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            label.append(int(dirs[5:]))
        labels.append(label)
        rg_0 = []
        pa_t = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)

            rg = get_rg()
            rg_0.append(rg)
            os.chdir("..")
        os.chdir("..")
        for dirs_p in p_dirs:
            rg_temp = []

            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)

                readInfo = open("Poly_Info.txt", "r")
                line = readInfo.readline()
                while "abs" not in line:
                    line = readInfo.readline()

                line_info = line.split()
                rg_temp.append(float(line_info[1]))

                #rg = get_rg()
                #rg_temp.append(rg)
                os.chdir("..")
            os.chdir("..")
            pa_t.append(np.array(rg_temp))
        pa_t = np.array(pa_t)
        pa_all.append(pa_t)
        os.chdir("..")
    sizes = []
    dist = parentFolders[0].split("_")[0]
    for dirs in parentFolders:
        sizes.append(int(dirs.split("_")[2]))


    fig, axs = plt.subplots(1,len(sizes))

    fig.suptitle("Percent Adsorption of "+dist)


    for i in range(len(labels)):
        labels[i] = np.around(np.array(labels[i])/labels[i][len(labels[i])-1],2)

    fig.set_size_inches(17, 8)
    index = 0
    pic = np.arange(0,1,.1)
    for array in pa_all:
        label_spaced = []
        j = 0
        for i in range(0,101):
            if (i/100== labels[index][j]):
                label_spaced.append(str(labels[index][j]))
                j = j +1
            else:
                label_spaced.append("")

        x = axs[index].imshow(normalizeHeatmaps(array,labels[index]),vmin=0,vmax=1)
        axs[index].set_xticks(np.arange(len(label_spaced)), labels=label_spaced)
        axs[index].set_yticks(np.arange(len(pic)), labels=np.around(pic,1))
        axs[index].title.set_text(str(sizes[index])+" Beads")
        axs[index].set_xlabel("Hinge Fraction")
        axs[index].set_ylabel("PiC")
        plt.colorbar(x,ax=axs[index],shrink=.5)
        axs[index].set_aspect(10)
        index = index+1
    fig.show()




def createRgHeatMap(parentFolders):
    rg_all = []
    labels = []
    for parentFolder in parentFolders:
        os.chdir(os.getcwd()+"\\"+parentFolder)
        p_dirs = next(os.walk('.'))[1]
        x_vals = []
        for dirs in p_dirs:
            x_vals.append(float(dirs[2:]))
        os.chdir(os.getcwd() + "\\" + p_dirs[0])
        c_dirs = next(os.walk('.'))[1]
        label = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            label.append(int(dirs[5:]))
        labels.append(label)
        rg_0 = []
        rg_t = []
        c_dirs.sort(key = lambda x: int(x.split("_")[2]))
        for dirs in c_dirs:
            os.chdir(os.getcwd() + "\\" + dirs)

            rg = get_rg()
            rg_0.append(rg)
            os.chdir("..")
        os.chdir("..")
        rg_0 = np.array((rg_0))

        for dirs_p in p_dirs:
            rg_temp = []

            os.chdir(os.getcwd() + "\\" + dirs_p)
            c_dirs = next(os.walk('.'))[1]
            c_dirs.sort(key = lambda x: int(x.split("_")[2]))
            for dirs in c_dirs:
                os.chdir(os.getcwd() + "\\" + dirs)
                rg = get_rg()
                rg_temp.append(rg)
                os.chdir("..")
            os.chdir("..")
            rg_t.append(np.array(rg_temp)/rg_0)
        rg_t = np.array(rg_t)
        rg_all.append(rg_t)
        os.chdir("..")

    sizes = []
    dist = parentFolders[0].split("_")[0]
    for dirs in parentFolders:
        sizes.append(int(dirs.split("_")[2]))


    fig, axs = plt.subplots(1,len(sizes))

    fig.suptitle("Radius of Gyration of "+dist)


    for i in range(len(labels)):
        labels[i] = np.around(np.array(labels[i])/labels[i][len(labels[i])-1],2)


    index = 0
    pic = np.arange(0,1,.1)
    fig.set_size_inches(17, 8)
    for array in rg_all:
        label_spaced = []
        j = 0
        for i in range(0,101):
            if (i/100== labels[index][j]):
                label_spaced.append(str(labels[index][j]))
                j = j +1
            else:
                label_spaced.append("")

        x = axs[index].imshow(normalizeHeatmaps(array,labels[index]),vmin=0,vmax=1.2)
        axs[index].set_xticks(np.arange(len(label_spaced)), labels=label_spaced)
        axs[index].set_yticks(np.arange(len(pic)), labels=np.around(pic,1))
        axs[index].title.set_text(str(sizes[index])+" Beads")
        axs[index].set_xlabel("Hinge Fraction")
        axs[index].set_ylabel("PiC")
        plt.colorbar(x,ax=axs[index],shrink=.5)
        axs[index].set_aspect(10)
        index = index+1
    fig.show()

