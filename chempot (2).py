#!/home/phys/qif/anaconda3/bin/python

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 07:01:11 2021

@author: dugue
"""

import sys
import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester

def matproj_formation_energy(formula):
    prop=['pretty_formula','material_id','energy_per_atom','formation_energy_per_atom']
    with MPRester("VUswQqBEWe4VFBZD25") as m:
        prop_dict={}
        for i in m.get_data(formula):
            if not i['pretty_formula'] in prop_dict or i['formation_energy_per_atom']<prop_dict[i['pretty_formula']][-1]:
                prop_dict[i['pretty_formula']]=[i[j] for j in prop[1:]]
    
    return prop_dict

def formula2dict(formula):
    chem_dict={}
    ele=""
    num=0
    for i in formula:
        if i.isupper():
            if ele:
                if not num:
                    chem_dict[ele]=1
                else:
                    chem_dict[ele]=num
            ele=i
            num=0
        elif i.islower():
            ele+=i
        elif i.isdigit:
            num=num*10+int(i)
        else:
            print(i)
    
    if not num:
        chem_dict[ele]=1
    else:
        chem_dict[ele]=num
        
    return chem_dict

abc=[[1/3**0.5,1],[0,0],[2/3**0.5,0]]
labc=[[0,0],[-3**0.5,2],[3**0.5,0]]

def parallel_line(ref_line,distance):
    temp=[ref_line[0],ref_line[1]+distance*(ref_line[0]**2+1)**0.5 if ref_line==labc[0] else ref_line[1]-distance*(ref_line[0]**2+1)**0.5]
    return temp

def crosspoint(line1,line2):
    if abs(line1[0]-line2[0])<1e-2:
        temp=[-(line1[1]-line2[1])/1e-2,(line1[0]*line2[1]-line2[0]*line1[1])/1e-2]
    else:
        temp=[-(line1[1]-line2[1])/(line1[0]-line2[0]),(line1[0]*line2[1]-line2[0]*line1[1])/(line1[0]-line2[0])]
    return temp

def distancetoline(line,xy):
    temp=(xy[0]*line[0]+line[1]-xy[1])/(line[0]**2+1)**0.5
    return abs(temp)

def tritoxy(tri):
    if type(tri[0])==bool:
        xy=crosspoint(parallel_line(labc[1],tri[1]),parallel_line(labc[2],tri[2]))
    elif type(tri[1])==bool:
        xy=crosspoint(parallel_line(labc[0],tri[0]),parallel_line(labc[2],tri[2]))
    else:
        xy=crosspoint(parallel_line(labc[0],tri[0]),parallel_line(labc[1],tri[1]))
    return xy

def xytotri(xy):
    return [round(distancetoline(i,xy),3) for i in labc]

formation_mat=[[2,1,6,9*-1.92235],[1,0,1,2*-2.25495],[0,1,4,5*-1.4272],[0,1,1,2*-0.81084],[0,1,2,3*-1.40226]]

def readtri(tri,mat):
    return [round(i*j/mat[-1],3) if not type(i)==bool else i for i,j in zip(tri,mat[:-1])]

def printtri(tri,mat):
    return [round(i/j*mat[-1],3) for i,j in zip(tri,mat[:-1])]

def line2point(point1,point2):
    if abs(point1[0]-point2[0])<1e-2:
        temp=[(point1[1]-point2[1])/1e-2,(point1[0]*point2[1]-point2[0]*point1[1])/1e-2]
    else:
        temp=[(point1[1]-point2[1])/(point1[0]-point2[0]),(point1[0]*point2[1]-point2[0]*point1[1])/(point1[0]-point2[0])]
    return temp

def plot(formula):
    compound0=matproj_formation_energy(formula)
    # {'Cs2SnCl6': ['mp-608555', -3.3898603444444446, -1.922349615651341]}
    compound={}
    # {'CsSn': ['mp-571056', -2.760116853125, -0.3076297395905172],
    #  'Cs4Sn23': ['mp-2496', -3.7147001966666666, -0.16648678358237556],
    #  'CsSn3': ['mp-865565', -3.38292671375, -0.15189722948275897],
    #  'CsCl': ['mp-573697', -3.31991964, -2.2549501239655174],
    #  'SnCl': ['mp-978928', -3.4328984175, -0.8108441600000003],
    #  'SnCl2': ['mp-569152', -3.561809344166667, -1.4022609525],
    #  'SnCl4': ['mp-29866', -3.2167423900000003, -1.4271986910000003]}
    dict0=formula2dict(formula)
    # {'Cs': 2, 'Sn': 1, 'Cl': 6}
    
    # rename formula
    for formula in compound0:
        pass
    
    compound_set=[]
    for i in dict0:
        for j in dict0:
            if not i==j and not {i,j} in compound_set:
                compound.update(matproj_formation_energy(i+"-"+j))
                compound_set.append({i,j})
        
    order=list(dict0.keys())
    mat0=[dict0[i] for i in order]
    mat0.append(sum(dict0.values())*compound0[formula][-1])
    
    for i in range(3):
        print(order[i]+" "+str(printtri(xytotri(abc[i]),mat0)))
        plt.annotate(order[i],xy=abc[i],xytext=(0,0),textcoords='offset points',fontsize=14)
    
    for i in range(3):
        for j in range(3)[i+1:]:
            temp=list(zip(abc[i],abc[j]))
            plt.plot(temp[0],temp[1],color="k")
            
    #temp=list(zip(*abc))
    #plt.scatter(temp[0],temp[1])
    lines={}
    for i in compound:
        dict1=formula2dict(i)
        # {'Cs': 1,'Cl': 1}
        tri=[]
        
        for j in dict1:
            tri.append([])
            for k in order:
                if k in dict1:
                    if j==k:
                        tri[-1].append(sum(dict1.values())*compound[i][-1]/dict1[k])
                    else:
                        tri[-1].append(0)
                else:
                    tri[-1].append(False)
        
        temp0=tritoxy(readtri(tri[0],mat0))
        temp1=tritoxy(readtri(tri[1],mat0))
        
        lines[i]=line2point(temp0,temp1)
        
        
        plt.plot([temp0[0],temp1[0]],[temp0[1],temp1[1]])
        
        #print(i,printtri(xytotri(temp0),mat0),printtri(xytotri(temp1),mat0))

        plt.annotate(i,xy=temp0,xytext=(0,0),textcoords='offset points',fontsize=14)
        #plt.annotate(i,xy=temp1,xytext=(0,0),textcoords='offset points',fontsize=14)
    
    for i,j in zip(order,labc):
        lines[i+"_rich"]=j
    
    line_set=[]
    for i in lines:
        for j in lines:
            if not i==j and not {i,j} in line_set:
                print(i,j,printtri(xytotri(crosspoint(lines[i],lines[j])),mat0))
                line_set.append({i,j})
    
    ax=plt.gca()
    ax.set_aspect(1)
    plt.savefig("pot.png",dpi=300)
    return None

if __name__ == '__main__':
    for i in sys.argv[1:]:
        plot(i)
