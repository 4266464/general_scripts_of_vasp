# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 03:14:34 2021

@author: dugue
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def plot_tdm(path, label_path, fig_path):
    with open(path, "r") as f:
        lines = f.readlines()

    temp_list=[]
    for i in lines[2:]:
        temp_list.append([float(j) for j in i.split(" ") if j][0:5:4][:2])
    temp_list=list(zip(*temp_list))

    #model = interp1d(temp_list[0],temp_list[1],kind="cubic")
    #x=np.linespace(temp_list[0][0],temp_list[0][-1],500)
    #y=model(x)
    plt.plot(temp_list[0],temp_list[1])

    with open(label_path, "r") as f:
        lines = f.readlines()

    label = []
    value = []
    for i in lines[1:]:
        if len(i) > 2:
            temp = [j for j in i.replace("\n", "").split(" ") if j]
            label.append(temp[0])
            value.append(float(temp[1]))
            plt.plot([float(temp[1]), float(temp[1])], [0, 30], "b--")


        else:
            break

    plt.xticks(value, label,fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlim(min(value), max(value))
    plt.ylim(0, 3)
    plt.ylabel("$\leftangle D \\rightangle^2$ (Debye$^2$)", fontsize=18)
    #plt.title(path)
    plt.savefig(fig_path + ".png", dpi=300)
    #plt.show()
    #plt.close('all')
    return None


# plot_tdm("TDM_DW_Zr.dat",
#           "KLABELS_Zr", "tdm_dw_Zr")
#plot_tdm("TDM_SOC_nb.dat","KLABELS_nb", "tdm_nb")
plot_tdm("hse06_TDM.dat","hse06_KLABELS", "hse06_tdm")