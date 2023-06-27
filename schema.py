# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 22:21:58 2021

@author: dugue
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist

def xy2(x,y,a=0.1,r=4,n=100):
    
    lx=np.linspace(x-r,x+r,n)
    ly=a*(lx-x)**2+y
    return [lx,ly]

grd=np.array([0,0])
sp=np.array([1,2])
cb=np.array([-4,3])
vb=np.array([4,2.5])

fgrd=lambda x:0.05*x**2
fcb=lambda x:0.1*(x--4)**2+3
f_cb=lambda x:0.1*(x--4)**2+3-0.2
fvb=lambda x:0.1*(x-4)**2+2.5
f_vb=lambda x:0.1*(x-4)**2+2.5-0.2
fsp=lambda x:0.1*(x-1)**2+2
f_sp=lambda x:0.1*(x-1)**2+2+0.2
    


def plott():
    
    # fig = plt.figure()
    # #使用axisartist.Subplot方法创建一个绘图区对象ax
    # ax = axisartist.Subplot(fig, 111)  
    # #将绘图区对象添加到画布中
    # fig.add_axes(ax)
    # #通过set_visible方法设置绘图区所有坐标轴隐藏
    # ax.lattice[:].set_visible(False)
    
    # #ax.new_floating_axis代表添加新的坐标轴
    # ax.lattice["x"] = ax.new_floating_axis(0,-1)
    # #给x坐标轴加上箭头
    # ax.lattice["x"].set_axisline_style("->", size = 2.0)
    # #添加y坐标轴，且加上箭头
    # ax.lattice["y"] = ax.new_floating_axis(1,-8)
    # ax.lattice["y"].set_axisline_style("->", size = 2.0)
    # #设置x、y轴上刻度显示方向
    # ax.lattice["x"].set_axis_direction("top")
    # ax.lattice["y"].set_axis_direction("right")
    # ax.lattice[:].major_ticks.set_tick_out(False)
    # ax.lattice["x"].label.set_text("configurations")
    # ax.lattice["y"].label.set_text("energy")
    for i in [sp,vb,cb]:
        plt.plot(*xy2(*i),"k")
    
    plt.plot(*xy2(*grd,a=0.05,r=5),"k")
    
    _cb=cb+[0,-0.2]
    _vb=vb+[0,-0.2]
    _sp=sp+[0,0.2]
    for i in [_sp,_cb,_vb]:
        plt.plot(*xy2(*i,r=1),"r--")
    
    uarrow={(0,fgrd(0)):[(0,fsp(0)),(0,fcb(0)),(0,fvb(0))]}
    darrow={(-4,f_cb(-4)):[(-4,fgrd(-4))],
            (4,f_vb(4)):[(4,fgrd(4))],
            (1,f_sp(1)):[(1,fgrd(1))]}
    
    for i in uarrow:
        for j in uarrow[i]:
            plt.plot([i[0],j[0]],[i[1],j[1]],"b-^")
    
    for i in darrow:
        for j in darrow[i]:
            plt.plot([i[0],j[0]],[i[1],j[1]],"r-v")
    
    double_arrow={(0,fgrd(0)):[(-4,fcb(-4)),(4,fvb(4)),(1,fsp(1))]}
    
    for i in double_arrow:
        for j in double_arrow[i]:
            plt.plot([i[0],j[0]],[i[1],j[1]],"c--")
    
    word_size=16
    plt.annotate("grd",xy=[5,1],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("cb",xy=[-7,4],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("$sp$",xy=[-2,3],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("vb",xy=[8,4],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("$^3P_1$",xy=[0.5,2.4],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    
    #coef=np.polyfit((1-1e-5,1+1e-5,2),(pos[1][0],pos[1][0],pos[1][1]),2)
    #plt.plot(x,np.poly1d(coef)(x))
    #plt.annotate(replace_dict[i],xy=energy_dict[i],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.xticks([])
    plt.yticks([])
    plt.plot([-10,10],[-0.5,-0.5],"k->")
    plt.plot([-9,-9],[-1,5],"k-^")
    
    plt.xlabel('Configurations',fontsize=18)
    #plt.ylim(grd-.6*(ex-grd),ex+.6*(ex-grd))
    plt.ylabel('Energy',fontsize=18)
    #plt.lattice("off")
    
    #ax.spines["left"].set_color("none")
    #ax.spines["bottom"].set_color("none")
    #plt.title("configuration coordinate curves")
    #plt.legend()
    #plt.show()
    plt.savefig("show1.png",dpi=300)
    plt.close('all')
    return None



def plott2():
    
    lines=[[[0.5,3.5],[0,0]],
           [[0.5,3.5],[4,4]],
           [[0.5,1.5],[3.5,3.5]],
           [[1,3.5],[0.5,0.5]],
           [[1.75,2.75],[3.3,3.3]],]
    
    arrow=[[[0.75,0.75],[0,3.5]],
           [[1.25,1.25],[0.5,3.5]],
           [[3.25,3.25],[0.5,4]],
           [[2.25,2.25],[0.5,3.3]]]
    
    
    plt.plot([-2.55,-2.55],[0.5,3.6],"b-^")
    plt.plot([-2.35,-2.75],[3.6,3.6],"b--")
    plt.plot([-1.95,-1.95],[0.8,3.3],"r-v")
    plt.plot([-1.75,-2.15],[0.8,0.8],"r--")
    word_size=16
    plt.annotate("VBM",xy=[-3.5,-0.35],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("CBM",xy=[-3.5,4.1],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("3+/4+",xy=[-3.2,0.6],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("3+*/4+",xy=[-2.3,3.45],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("2+/3+",xy=[-1.3,3.6],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("MMCT",xy=[-3.2,1.8],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    #plt.annotate("ex sp",xy=[-2.76,1.2],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("$sp$",xy=[-1.9,1.8],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    #plt.annotate("em sp",xy=[-2.2,2.4],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("IV",xy=[-1.2,1.8],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("CT",xy=[-0.7,1.8],xytext=(0,0),textcoords='offset points',fontsize=word_size)
    
    for i in lines:
        i[0]=-np.array(i[0])
        plt.plot(*i,"k")
    
    for i in arrow:
        i[0]=-np.array(i[0])
        plt.plot(*i,"c--^")
    #coef=np.polyfit((1-1e-5,1+1e-5,2),(pos[1][0],pos[1][0],pos[1][1]),2)
    #plt.plot(x,np.poly1d(coef)(x))
    #plt.annotate(replace_dict[i],xy=energy_dict[i],xytext=(0,0),textcoords='offset points',fontsize=word_size)
          
    #plt.xlabel('Configurations')
    #plt.ylim(grd-.6*(ex-grd),ex+.6*(ex-grd))
    #plt.ylabel('Energy')
    plt.ylim(-0.4,4.5)
    plt.xticks([])
    plt.yticks([])
    #plt.title("charge transition level")
    # plt.legend()
    plt.axis("off")
    #plt.show()
    plt.savefig("show2.png",dpi=300)
    plt.close('all')
    return None

plott()
#plott2()