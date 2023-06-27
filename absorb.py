# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 17:06:23 2021

@author: dugue
"""
import matplotlib.pyplot as plt

def absorb(path):
    xy=[]
    with open(path,"r") as f:
        for i in f.readlines():
            if "#" in i:
                xy.append([j for j in i.replace("#","").split(" ") if j][:2])
            else:
                xy.append([float(j) for j in i.split(" ") if j][:2])
            if len(xy)>1 and xy[-1][0]>10:
                break
    
    for i in range(len(xy))[3:-1]:
        if (xy[i][1]-xy[i-1][1])>0 and (xy[i+1][1]-xy[i][1])<0:
            plt.scatter(*xy[i],color="red")
            plt.annotate(str(xy[i][0]),xy[i])
        # if (xy[i][1]-xy[i-1][1])-(xy[i-1][1]-xy[i-2][1])<0 and (xy[i+1][1]-xy[i][1])-(xy[i][1]-xy[i-1][1])>0:
        #     plt.scatter(*xy[i],color="blue")
        #     plt.annotate(str(xy[i][0]),xy[i])
        
    xy=list(zip(*xy))
    plt.plot(xy[0][1:],xy[1][1:])
    plt.xlabel(xy[0][0][0].upper()+xy[0][0][1:]+" (eV)")
    plt.ylabel(xy[1][0])
    plt.xlim(0,)
    plt.ylim(0,)
    plt.show()
    #plt.savefig(path.split(".")[0]+".png",dpi=300)
    #plt.close('all')
    return None

#absorb("pbe_new_ABSORB.dat")
absorb("hse06_old_ABSORB.dat")
#absorb("hse06_new_ABSORB.dat")
    