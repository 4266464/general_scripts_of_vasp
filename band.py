# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 03:14:34 2021

@author: dugue
"""

"""
readme

three type band structure file supported accompany with k path label file
the band structure file looks like:
1. classic
#K-Path     Energy-Level
#Band-index    1
  0.000      -41.1290
  0.024      -41.1290
2. spin-up/spin-down
#K-Path             Spin-Up       Spin-Down
#Band-index    1
  0.000      -36.1400      -36.1400
  0.024      -36.1400      -36.1400
3. fat band
#K-Path    Energy     s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
#Band-index    1
   0.000   -23.679  0.000  0.000  0.000  0.000  0.195  0.195  0.199  0.195  0.199  0.984
   0.030   -23.679  0.000  0.000  0.000  0.000  0.183  0.183  0.211  0.219  0.187  0.984

the k path label file looks like
K-Label    K-Coordinate in band-structure plots 
GAMMA              0.000
X                  0.459
"""
import matplotlib.pyplot as plt
import numpy as np
def plot_band(band_path, label_path, fig_path, select_band=False,energy_shift=0.):
    if select_band:
        with open(band_path[0][0] if type(band_path)==list else band_path, "r") as f:
            first_lines=f.readlines()

        print("band energy")
        flag = False
        for i in first_lines[1:]:
            if "#" in i:
                print(i.split(" ")[-1].replace("\n", ""), end=" ")
                flag = True
            elif flag:
                print([j for j in i.split(" ") if j][1])
                flag = False

        band_range=input("band range:")
        band_range=[float(i) for i in band_range.split(" ")]
    else:
        band_range = [-1, 10000]

    band_max=[]

    # multi bandstructure file / fat band
    if type(band_path) == list:
        lines = []
        for i in band_path:
            with open(i[0], "r") as f:
                lines.append(f.readlines())
        for k, l in zip(lines, band_path):
            keys = [i.replace("x2-y2", "dx2-y2") for i in k[0].replace("#", "").replace("\n", "").split(" ") if i]
            temp = []
            x = []
            y = []
            flag = False
            band_count=0
            for i in k[1:]:
                if "#" in i:
                    band_count+=1
                    if float(i.split(" ")[-1].replace("\n", "")) < band_range[0]:
                        continue
                    elif float(i.split(" ")[-1].replace("\n", "")) > band_range[1]:
                        break
                    else:
                        if band_max and y:
                            index_min = y.index(min(y))
                            if y[index_min]>band_max[1] and -0.5<band_max[1]<0.5:
                                #plt.plot([band_max[0], x[index_min]], [band_max[1]-energy_shift, y[index_min]-energy_shift], linestyle="dashdot",marker="o")
                                print(band_count, band_max[1], y[index_min],y[index_min]-band_max[1])
                        if y:
                            index_max=y.index(max(y))
                            band_max=[x[index_max],y[index_max]]
                        flag = True
                        temp = []
                        x = []
                        y = []
                elif flag:
                    curr = [float(j) for j in i.replace("#", "").replace("\n", "").split(" ") if j]
                    if len(curr)==0 or (len(x) > 1 and curr[0] == x[-1]):
                        pass
                    # 剔除不连续点
                    elif len(x) > 2 and abs(
                            (curr[1] - y[-1]) / (curr[0] - x[-1]) - (y[-1] - y[-2]) / (x[-1] - x[-2])) > 5:
                        # plt.plot(x[-1],y[-1],"rx")
                        x[-1] = curr[0]
                        y[-1] = curr[1]
                    else:
                        x.append(curr[0])
                        y.append(curr[1])
                        curr = [curr[0], curr[1], sum([curr[j] for j in range(len(curr)) if keys[j][0] == l[1]])]
                        if temp:
                            plt.plot([temp[0], curr[0]], [temp[1]-energy_shift, curr[1]-energy_shift], linewidth=(temp[2] + curr[2]), color=l[2])
                        temp = curr[:]

    # classic and spin-up/down
    else:
        with open(band_path, "r") as f:
            lines=f.readlines()
        x = []
        y = []
        flag = False
        for i in lines[1:]:
            if "#" in i:
                if float(i.split(" ")[-1].replace("\n", "")) < band_range[0]:
                    continue
                elif float(i.split(" ")[-1].replace("\n", "")) > band_range[1]+1:
                    # if x:
                    #     if type(y[0])==list:
                    #         y1,y2=list(zip(*y))
                    #         plt.plot(x, y1, color="r")
                    #         plt.plot(x, y2, color="b")
                    #     else:
                    #         y= list(zip(*y))
                    #         plt.plot(x, y, color="b")
                    break
                else:
                    flag = True
                    if x:
                        if type(y[0]) == list:
                            y1, y2 = list(zip(*y))
                            #plt.plot(x, y1, color="r",linewidth=2)
                            plt.plot(x, np.array(y2)-energy_shift if min(np.array(y2))>0 else np.array(y2), color="b")
                        else:
                            #y = list(zip(*y))
                            plt.plot(x, np.array(y)-energy_shift if min(np.array(y))>0 else np.array(y), color="b")
                        x = []
                        y = []
            elif flag and len(i) > 5:
                curr = [float(j) for j in i.replace("#", "").replace("\n", "").split(" ") if j]
                if len(x) > 1 and curr[0] == x[-1]:
                    pass
                # check non-continuous points
                # elif len(x) > 2 and abs((curr[1] - y[-1]) / (curr[0] - x[-1]) - (y[-1] - y[-2]) / (x[-1] - x[-2])) > 10:
                #     # plt.plot(x[-1],y[-1],"rx")
                #     x[-1] = curr[0]
                #     y[-1] = curr[1:]
                else:
                    x.append(curr[0])
                    y.append(curr[1:] if len(curr)>2 else curr[1])

    with open(label_path, "r") as f:
        lines = f.readlines()
    label = []
    value = []
    for i in lines[1:]:
        if len(i) > 2:
            temp = [j for j in i.replace("\n", "").split(" ") if j]
            #label.append(temp[0].replace("GAMMA","$\Gamma$"))
            label.append(temp[0])
            value.append(float(temp[1]))
            plt.plot([float(temp[1]), float(temp[1])], [-40, 40], "b--")

        else:
            break

    plt.xticks(value, label,fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(min(value), max(value))
    plt.ylim(-4, 10)
    #plt.xlabel(x[0])
    plt.ylabel("Energy (eV)",fontsize=22)
    # plt.title(band_path)
    plt.savefig(fig_path + ".png", dpi=300)
    #plt.show()
    plt.close('all')
    return None

if __name__ == '__main__':
    # plot_band([["hse06_PBAND_Zr.dat", "d", "green"],
    #            ["hse06_PBAND_Cl.dat", "p", "purple"]],
    #           "hse06_KLABELS", "hse06_band")

    # plot_band([["PBAND_Zr_UP_Zr.dat", "d", "green"],
    #            ["PBAND_Cl_UP_Zr.dat", "p", "purple"]],
    #           "KLABELS_Zr", "pbe_band",energy_shift=-0.265)
    #plot_band("7BAND.dat","7KLABELS", "cao")
    # plot_band("BAND.dat_Sn", "KLABELS_Sn", "Cs2SnCl6")
    # plot_band("BAND.dat_Zr", "KLABELS_Zr", "Cs2ZrCl6")
    # plot_band("BAND.dat_Hf", "KLABELS_Hf", "Cs2HfCl6")
    # plot_band("BAND.dat_Sn_soc", "KLABELS_Sn", "Cs2SnCl6_soc")
    # plot_band("BAND.dat_Zr_soc", "KLABELS_Zr", "Cs2ZrCl6_soc")
    #plot_band("BAND.dat_Zr", "KLABELS_Zr", "Cs2ZrCl6",energy_shift=-1)
    plot_band("BAND.dat", "KLABELS2", "xhm",energy_shift=0)
