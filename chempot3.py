# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 02:14:04 2021

@author: dugue
"""

import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
import numpy as np
import json
import random
import os

def matproj_formation_energy(formula):
    prop=['unit_cell_formula','material_id','energy_per_atom','formation_energy_per_atom']
    get_data=False
    #flag=False
    with open("mpb","r") as f:
        for i in f.readlines():
            if type(get_data)==list:
                if "search " in i:
                    break
                else:
                    get_data.append(json.loads(i.replace("\n","")))
            elif i.replace("search ","").replace("\n","")==formula:
                get_data=[]
                #flag=True
    
    if type(get_data)==bool:
        get_data=[]
        with open("mpb","a") as f:
            with MPRester("VUswQqBEWe4VFBZD25") as m:
                #temp=m.get_data(formula)
                #if temp:
                f.write("search "+formula+"\n")
                for i in m.get_data(formula):
                    get_data.append(i)
                    f.write(json.dumps(i)+"\n")
                
                
    prop_dict={}
    for i in get_data:
        agcd=gcd(list(i['unit_cell_formula'].values()))
        if agcd:
            for j in i['unit_cell_formula']:
                i['unit_cell_formula'][j]/=agcd
        item=tuple(sorted(list(i['unit_cell_formula'].items()),key=lambda x:x[0]))
        if not item in prop_dict or i['formation_energy_per_atom']<prop_dict[item][-1]:
            prop_dict[item]=[i[j] for j in prop[1:]]
    
    return prop_dict

def gcd(alist):
    blist=list(set(alist[:]))
    clist=blist[:]
    while True:
        if 1 in blist:
            return False
        else:
            for i in blist:
                for j in blist:
                    sub=abs(i-j)
                    if sub and not sub in clist:
                        clist.append(sub)
            if len(clist)==len(blist):
                return min(blist)
            else:
                blist=clist[:]

def formula2tuple(formula):
    chem_dict={}
    ele=""
    num=0.
    for i in formula:
        if i.isupper():
            if ele:
                if not num:
                    chem_dict[ele]=1.
                else:
                    chem_dict[ele]=num
            ele=i
            num=0.
        elif i.islower():
            ele+=i
        elif i.isdigit:
            num=num*10+float(i)
        else:
            print("illegal "+i)
    
    if not num:
        chem_dict[ele]=1.
    else:
        chem_dict[ele]=num
    chem_dict=tuple(sorted(list(chem_dict.items()),key=lambda x:x[0]))
    
    return chem_dict

def tuple2formula(chem_tuple):
    formula=""
    for i in chem_tuple:
        for j in i:
            if type(j)==float:
                if j==1:
                    pass
                else:
                    formula+=str(int(j))
            else:
                formula+=j
    return formula

def line2point(point1,point2):
    if abs(point1[0]-point2[0])<1e-5:
        temp=[(point1[1]-point2[1])/1e-5,(point1[0]*point2[1]-point2[0]*point1[1])/1e-5]
    else:
        temp=[(point1[1]-point2[1])/(point1[0]-point2[0]),(point1[0]*point2[1]-point2[0]*point1[1])/(point1[0]-point2[0])]
    return temp

# abc=[[1/3**0.5,1],[0,0],[2/3**0.5,0]]
# labc=[line2point(*i) for i in [[j for j in abc if not j==i] for i in abc]]

# def parallel_line(ref_line,distance):
#     if ref_line==labc[0]:
#         return [ref_line[0],ref_line[1]+distance*(ref_line[0]**2+1)**0.5]
#     else:
#         return [ref_line[0],ref_line[1]-distance*(ref_line[0]**2+1)**0.5]

# def distancetoline(line,xy):
#     temp=(xy[0]*line[0]+line[1]-xy[1])/(line[0]**2+1)**0.5
#     return abs(temp)


#abc=[[1,1],[0,1],[1,0]]
abc=[[0,0],[-1,0],[0,-1]]
labc=[line2point(*i) for i in [[j for j in abc if not j==i] for i in abc]]

def parallel_line(ref_line,distance):
    if ref_line==labc[0]:
        return [ref_line[0],ref_line[1]+distance*(ref_line[0]**2+1)**0.5/2**0.5]
    elif ref_line==labc[1]:
        return [ref_line[0],ref_line[1]+distance*(ref_line[0]**2+1)**0.5]
    else:
        return [ref_line[0],ref_line[1]-distance*(ref_line[0]**2+1)**0.5]

def distancetoline(line,xy):
    if line[0]==1:
        temp=(xy[0]*line[0]+line[1]-xy[1])/(line[0]**2+1)**0.5*2**0.5
    else:
        temp=(xy[0]*line[0]+line[1]-xy[1])/(line[0]**2+1)**0.5
    return abs(temp)

def crosspoint(line1,line2):
    if abs(line1[0]-line2[0])<1e-5:
        temp=[-(line1[1]-line2[1])/1e-5,(line1[0]*line2[1]-line2[0]*line1[1])/1e-5]
    else:
        temp=[-(line1[1]-line2[1])/(line1[0]-line2[0]),(line1[0]*line2[1]-line2[0]*line1[1])/(line1[0]-line2[0])]
    return temp

def tritoxy(tri):
    if type(tri[0])==bool:
        xy=crosspoint(parallel_line(labc[1],tri[1]),parallel_line(labc[2],tri[2]))
    elif type(tri[1])==bool:
        xy=crosspoint(parallel_line(labc[0],tri[0]),parallel_line(labc[2],tri[2]))
    else:
        xy=crosspoint(parallel_line(labc[0],tri[0]),parallel_line(labc[1],tri[1]))
    return xy

def xytotri(xy):
    temp=np.array([round(distancetoline(i,xy),3) for i in labc])
    temp/=np.linalg.norm(temp,ord=2)
    return temp

def readtri(tri,formula_tuple,formation_energy):
    return [i*j[1]/formation_energy if not type(i)==bool else i for i,j in zip(tri,formula_tuple)]

def printtri(tri,formula_tuple,formation_energy):
    return [round(i/j[1]*formation_energy,3) for i,j in zip(tri,formula_tuple)]



def condition(chempot_list,formation_list):
    return sum([i*j for i,j in zip(chempot_list,formation_list)])<formation_list[-1]

def chempot_sample(formula,fixed={},start_point=[],other_element=[],small=False):
    if False:
        pass
    # plot from csv
    # if os.path.isfile(formula+"_".join([i+"_"+str(round(fixed[i],3)) for i in fixed])+".csv"):
    #     with open(formula+"_".join([i+"_"+str(round(fixed[i],3)) for i in fixed])+".csv","r") as f:
    #         flag=False
    #         for i in f.readlines():
    #             if "sample" in i and not flag:
    #                 flag=True
    #             elif "range" in i:
    #                 break
    #             elif flag:
    #                 plt.scatter([float(j) for j in i.split(",")])
                
                    
            
    else:
        with open(formula+"_".join([i+"_"+str(round(fixed[i],3)) for i in fixed])+".csv","w") as f:
            # "Cs2SnCl6"
            formula_tuple=formula2tuple(formula)
            # SOC crystal and mocular
            # soc={"Bi":[-0.473,-0.808],
            #      "Cs":[-0.131,-0.143],
            #      "Na":[-0.000,-0.000],
            #      "Cl":[-0.007,-0.002]}
            soc={}
    
            # (('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0))
            formation_mat={}
            other_formation_mat={}
            for i in range(2**(len(formula_tuple)+len(other_element)))[1:]:
                code=bin(i)
                code=code.replace("0b","".join(["0" for j in range(len(formula_tuple)+len(other_element)-len(code)+2)]))
                # 011
                ions=[]
                for j,k in zip(code,list(dict(formula_tuple).keys())+other_element):
                    if j=="1":
                        ions.append(k)
                formation_dict=matproj_formation_energy("-".join(ions))
                # {(('Cl', 3.0), ('Cs', 1.0), ('Sn', 1.0)): ['mp-27394', -3.533219568, -1.8115027265862065], 
                # (('Cl', 5.0), ('Cs', 1.0), ('Sn', 2.0)): ['mp-30164', -3.551530498125, -1.6656268253663797],
                # (('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0)): ['mp-608555', -3.3898603444444446, -1.922349615651341]}
                for j in formation_dict:
                    temp_formation_mat=[]
                    formation_energy=sum(list(zip(*j))[1])*formation_dict[j][-1]
                    if soc:
                        for k in dict(j):
                            if k in soc:
                                formation_energy+=dict(j)[k]*(soc[k][0]-soc[k][1])
                    
                    for k in list(dict(formula_tuple).keys())+other_element:
                        if k in dict(j):
                            if k in fixed:
                                formation_energy-=dict(j)[k]*fixed[k]
                            else:
                                temp_formation_mat.append(dict(j)[k])
                        else:
                            if k in fixed or k in other_element:
                                pass
                            else:
                                temp_formation_mat.append(0)
                    
                    temp_formation_mat.append(formation_energy)
                    temp_formation_mat.append(tuple2formula(j))
                    
                    if set(other_element).intersection(set(dict(j))):
                        other_formation_mat[j]=temp_formation_mat[:]
                    else:
                        formation_mat[j]=temp_formation_mat[:]
            
            for i in other_formation_mat:
                print(i,other_formation_mat[i])
            print("----------------------------------------")
            for i in formation_mat:
                print(i,formation_mat[i])
            
            real_formula_tuple=tuple([i for i in formula_tuple if not i[0] in fixed])
            
            chempot_min=[formation_mat[formula_tuple][-2]/i for i in formation_mat[formula_tuple][:-2]]
            
            f.write(",".join([str(i[0]) for i in real_formula_tuple])+"\n")
            f.write(",".join([str(i) for i in chempot_min])+"\n")
                
            for i in formation_mat:
                f.write(",".join([str(j) for j in formation_mat[i]])+"\n")
            
            if len(real_formula_tuple)==3:
                line_mat={}
                ax1=plt.subplot(121,aspect=1)
                ax2=plt.subplot(122)
                for i in formation_mat:
                    if not i==formula_tuple and formation_mat[i][-2]<0:
                        num_notzero=len([1 for j in formation_mat[i][:-2] if not j==0])
                        if num_notzero==1:
                            line_mat[i]=[]
                            for j in range(3):
                                if not formation_mat[i][j]==0:
                                    line_mat[i].append(0)
                                    ratio=formation_mat[formula_tuple][j]/formation_mat[i][j]
                                else:
                                    line_mat[i].append(formation_mat[formula_tuple][j])
                            line_mat[i].append(formation_mat[formula_tuple][3]-ratio*formation_mat[i][3])
                            line_mat[i].append(formation_mat[i][4])
                        elif num_notzero==2:
                            line_mat[i]=formation_mat[i][:]
                        elif num_notzero==3:
                            for j in range(3):
                                ratio=formation_mat[formula_tuple][j]/formation_mat[i][j]
                                line_mat[i]=[]
                                for k in range(3):
                                    line_mat[i].append(formation_mat[formula_tuple][k]-ratio*formation_mat[i][k])
                                if len([1 for j in line_mat[i] if not j==0])==2:
                                    break
                            line_mat[i].append(formation_mat[formula_tuple][3]-ratio*formation_mat[i][3])
                            line_mat[i].append(formation_mat[i][4])
                
                lines={}
                for i in line_mat:
                    end_points=[]
                    for j in range(3):
                        if not line_mat[i][j]==0:
                            end_points.append([])
                            for k in range(3):
                                if j==k:
                                    end_points[-1].append(line_mat[i][3]/line_mat[i][k])
                                elif line_mat[i][k]==0:
                                    end_points[-1].append(False)
                                else:
                                    end_points[-1].append(0)
                    
                    end_points[0]=tritoxy(readtri(end_points[0],real_formula_tuple,formation_mat[formula_tuple][-2]))
                    end_points[1]=tritoxy(readtri(end_points[1],real_formula_tuple,formation_mat[formula_tuple][-2]))
                    lines[line_mat[i][4]]=line2point(end_points[0],end_points[1])
                    #ax1.plot(*list(zip(*end_points)),label=line_mat[i][4],linewidth=0.2)
                    ax1.plot(*list(zip(*end_points)),label=line_mat[i][4])
                    ax2.plot(*list(zip(*end_points)),label=line_mat[i][4])
                    #ax1.annotate(line_mat[i][4],xy=end_points[0],xytext=(0,0),textcoords='offset points',fontsize=6)
                    #ax1.annotate(line_mat[i][4],xy=end_points[1],xytext=(0,0),textcoords='offset points',fontsize=6)         
    
                for i,j in zip(formula_tuple,labc):
                    lines[i[0]+"_rich"]=j
                
                lines_set=[]
                for i in lines:
                    for j in lines:
                        if not i==j and not {i,j} in lines_set:
                            ##print(i,j,crosspoint(lines[i],lines[j]),printtri(xytotri(crosspoint(lines[i],lines[j])),real_formula_tuple,formation_mat[formula_tuple][-2]))
                            lines_set.append({i,j})
    
            if start_point:
                chempot_list=start_point
            else:
                for i in range(100000000):
                    res_formation=formation_mat[formula_tuple][-2]
                    chempot_list=[]
                    for j in real_formula_tuple:
                        ratio=random.random()
                        chempot_list.append(res_formation*ratio/j[1])
                        res_formation*=(1-ratio)
                    
                    for j in formation_mat:
                        if not j==formula_tuple:
                            if not condition(chempot_list,formation_mat[j][:-1]):
                                chempot_list=[]
                                break
                        
                    if chempot_list:
                        break
                
            if chempot_list:
                ax1.scatter(*tritoxy(readtri(chempot_list,real_formula_tuple,formation_mat[formula_tuple][-2])),color="red")
                #f.write("sample\n"+",".join([str(i) for i in chempot_list])+"\n")
                
                chempot_range=[[[],[]] for i in formation_mat[formula_tuple][:-2]]
                for j in range(len(chempot_list)):
                    chempot_range[j][0]=chempot_list
                    chempot_range[j][1]=chempot_list
                
                ratio_list=[i/j for i,j in zip(chempot_list,chempot_min)]
                count=0
                while count<100:
                    temp_ratio_list=[j+0.1*(2*random.random()-1) for j in ratio_list[:-1]]
                    temp_ratio_list.append(1-sum(temp_ratio_list))
                    for j in temp_ratio_list:
                        if j>1 or j<0:
                            temp_ratio_list=[]
                            break
                        else:
                            chempot_list=[i*j for i,j in zip(temp_ratio_list,chempot_min)]
                    
                    if not temp_ratio_list:
                        continue
                    
                    for j in formation_mat:
                        if not j==formula_tuple:
                            if not condition(chempot_list,formation_mat[j][:-1]):
                                chempot_list=[]
                                break
                        
                    if chempot_list:
                        count+=1
                        ratio_list=temp_ratio_list
                        #f.write(",".join([str(j) for j in chempot_list])+"\n")
                        
                        for j in range(len(chempot_list)):
                            if chempot_list[j]>chempot_range[j][0][j]:
                                chempot_range[j][0]=chempot_list
                                #print(real_formula_tuple[j][0],chempot_range[j])
                                
                            if chempot_list[j]<chempot_range[j][1][j]:
                                chempot_range[j][1]=chempot_list
                                #print(real_formula_tuple[j][0],chempot_range[j])
                        
                        
                        #if len(real_formula_tuple)==3:
                            #xy=tritoxy(readtri(chempot_list,real_formula_tuple,formation_mat[formula_tuple][-2]))
                            #ax1.scatter(*xy,color="red")
                            #ax2.scatter(*xy,color="red",s=0.1)
                
                f.write("range\n")
                
                res_dict={}
                for i,j in zip(real_formula_tuple,chempot_range):
                    for k,l in zip(j,["max","min"]):
                        if tuple(k) in res_dict:
                            res_dict[tuple(k)].append(i[0]+"_"+l)
                        else:
                            res_dict[tuple(k)]=[i[0]+"_"+l]
                print(res_dict)
                if other_element:
                    for i in res_dict:
                        other_min=[]
                        for j in other_formation_mat.values():
                            other_min.append([j[-1],(j[-2]-sum([i[k]*j[k] for k in range(len(i))]))/j[-3]])
                        other_min.sort(key=lambda x:x[-1])
                        #res_dict[i].append(other_min[0][0])
                        #res_dict[i].append(list(i)+[other_min[0][1]])
                        for j in other_min[:5]:
                            print(j[0],j[1])
                            f.write(j[0]+","+str(j[1])+"\n")
                    # else:
                    #     res_dict[i].append(list(i))
                    #     res_dict[i].append(tritoxy(readtri(i,real_formula_tuple,formation_mat[formula_tuple][-2])))
                xy_list=[]
                for i in res_dict:
                    xy_list.append([i,tritoxy(readtri(i,real_formula_tuple,formation_mat[formula_tuple][-2]))])
                    ax2.scatter(*xy_list[-1][1],color="blue")
                    ax2.annotate(",".join(res_dict[i]),xy_list[-1][1],xytext=(0,0),textcoords='offset points',fontsize=8)
                    
                    
                    # for j in res_dict[i][:-2]:
                    #     print(j,end=" ")
                    #     f.write(str(j)+",")
                    # print()
                    # f.write("\n")
                    # print(*res_dict[i][-2])
                    # f.write(",".join([str(j) for j in res_dict[i][-2]])+"\n")
                    

            
                if len(real_formula_tuple)==3:
                    lines=[]
                    for j,k in zip(real_formula_tuple,abc):
                        #print(k,j[0],printtri(xytotri(k),formula_tuple,formation_mat[0][-2]))
                        ax1.annotate(j[0]+"("+",".join([str(l) for l in printtri(xytotri(k),real_formula_tuple,formation_mat[formula_tuple][-2])])+")",
                                              k,xytext=(0,0),textcoords='offset points',fontsize=8)
                        # ax2.annotate(j[0]+str(printtri(xytotri(k),real_formula_tuple,formation_mat[formula_tuple][-2])),
                        #                       k,xytext=(0,0),textcoords='offset points',fontsize=4)
                        lines.append(k)
                   
                    for j in range(3):
                        for k in range(3)[j+1:]:
                            temp=list(zip(lines[j],lines[k]))
                            ax1.plot(temp[0],temp[1],color="k")
                            #ax1.plot(temp[0],temp[1],color="k",linewidth=0.2)
                            ax2.plot(temp[0],temp[1],color="k")
                            
                    #.lattice("off")
                    #ax=plt.gca()
                    #ax.set_aspect(1)
                    
                    #xy_list=[[i,res_dict[i][-1]] for i in res_dict]
                    xy_list.sort(key=lambda x:x[-1][0])
                    xmin=xy_list[0]
                    xmax=xy_list[-1]
                    xy_list.sort(key=lambda x:x[-1][1])
                    ymin=xy_list[0]
                    ymax=xy_list[-1]
                    
                    # xmax=[i[-1][0] for i in res_dict.values()]
                    # xmin=min(xmax)
                    # xmax=max(xmax)
                    # ymax=[i[-1][1] for i in res_dict.values()]
                    # ymin=min(ymax)
                    # ymax=max(ymax)
                    
                    
                    ax2.set_xlim(xmin[-1][0]-0.2*(xmax[-1][0]-xmin[-1][0]),xmax[-1][0]+0.2*(xmax[-1][0]-xmin[-1][0]))
                    ax2.set_ylim(ymin[-1][1]-0.2*(ymax[-1][1]-ymin[-1][1]),ymax[-1][1]+0.2*(ymax[-1][1]-ymin[-1][1]))
                    
                    ax1.set_xticks(np.linspace(-1,0,4))
                    ax1.set_xticklabels(np.round(np.linspace(printtri(xytotri([-1,0]),real_formula_tuple,formation_mat[formula_tuple][-2])[1],0,4),3))
                    ax1.set_yticks(np.linspace(0,-1,4))
                    ax1.set_yticklabels(np.round(np.linspace(0,printtri(xytotri([0,-1]),real_formula_tuple,formation_mat[formula_tuple][-2])[2],4),3))
                    ax1.tick_params(gridOn=True,grid_linewidth=0.5,
                                    bottom=False,top=True,left=False,right=True,
                                    labelbottom=False,labeltop=True,labelleft=False,labelright=True)
                    
                    
                    box=ax1.get_position()
                    ax1.set_position([box.x0,box.y0,box.width,box.height*0.8])
                    ax1.legend(fontsize=6,loc="lower left",bbox_to_anchor=[0,1.12],ncol=3)
                    
                    ax2.set_xticks(np.linspace(xmin[-1][0],xmax[-1][0],3))
                    ax2.set_xticklabels(np.round(np.linspace(xmin[0][1],
                                                    xmax[0][1],3),3))
                    ax2.set_yticks(np.linspace(ymin[-1][1],ymax[-1][1],4))
                    ax2.set_yticklabels(np.round(np.linspace(ymin[0][2],
                                                    ymax[0][2],4),3))
                    ax2.tick_params(gridOn=True,grid_linewidth=0.5,
                                    bottom=False,top=True,left=False,right=True,
                                    labelbottom=False,labeltop=True,labelleft=False,labelright=True)
                    #ax2.legend(fontsize=6,loc="upper right")
                    plt.savefig(formula+"_".join([i+"_"+str(round(fixed[i],3)) for i in fixed])+".png",dpi=300)
                    plt.close('all')
                    #plt.show()
            
            else:
                print("sample failed")
        return res_dict

chempot_sample("Cs2SnCl6")

# sample=False
# with open("Cs2NaBiCl6.csv","r") as f:
#     for i in f.readlines():
#         if "range" in i:
#             sample.sort(key=lambda x:x[1])
#             mid=len(sample)//2
#             for j in [0,mid,-1]:
#             #for j in [0]:
#                 chempot_sample("Cs2NaBiCl6",
#                                 fixed={"Cl":sample[j][1]},
#                                 start_point=[sample[j][0],sample[j][2],sample[j][3]],
#                                 other_element=["O"])
#             break
                
                
#         if type(sample)==list:
#             sample.append([float(j) for j in i.replace("\n","").split(",")])
            
#         if "sample" in i:
#             sample=[]
