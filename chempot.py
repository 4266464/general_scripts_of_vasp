# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 07:01:11 2021

@author: dugue
"""

import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
import numpy as np
import json

cnames={'aliceblue':'#F0F8FF',
'antiquewhite':'#FAEBD7',
'aqua':'#00FFFF',
'aquamarine':'#7FFFD4',
'azure':'#F0FFFF',
'beige':'#F5F5DC',
'bisque':'#FFE4C4',
'black':'#000000',
'blanchedalmond':'#FFEBCD',
'blue':'#0000FF',
'blueviolet':'#8A2BE2',
'brown':'#A52A2A',
'burlywood':'#DEB887',
'cadetblue':'#5F9EA0',
'chartreuse':'#7FFF00',
'chocolate':'#D2691E',
'coral':'#FF7F50',
'cornflowerblue':'#6495ED',
'cornsilk':'#FFF8DC',
'crimson':'#DC143C',
'cyan':'#00FFFF',
'darkblue':'#00008B',
'darkcyan':'#008B8B',
'darkgoldenrod':'#B8860B',
'darkgray':'#A9A9A9',
'darkgreen':'#006400',
'darkkhaki':'#BDB76B',
'darkmagenta':'#8B008B',
'darkolivegreen':'#556B2F',
'darkorange':'#FF8C00',
'darkorchid':'#9932CC',
'darkred':'#8B0000',
'darksalmon':'#E9967A',
'darkseagreen':'#8FBC8F',
'darkslateblue':'#483D8B',
'darkslategray':'#2F4F4F',
'darkturquoise':'#00CED1',
'darkviolet':'#9400D3',
'deeppink':'#FF1493',
'deepskyblue':'#00BFFF',
'dimgray':'#696969',
'dodgerblue':'#1E90FF',
'firebrick':'#B22222',
'floralwhite':'#FFFAF0',
'forestgreen':'#228B22',
'fuchsia':'#FF00FF',
'gainsboro':'#DCDCDC',
'ghostwhite':'#F8F8FF',
'gold':'#FFD700',
'goldenrod':'#DAA520',
'gray':'#808080',
'green':'#008000',
'greenyellow':'#ADFF2F',
'honeydew':'#F0FFF0',
'hotpink':'#FF69B4',
'indianred':'#CD5C5C',
'indigo':'#4B0082',
'ivory':'#FFFFF0',
'khaki':'#F0E68C',
'lavender':'#E6E6FA',
'lavenderblush':'#FFF0F5',
'lawngreen':'#7CFC00',
'lemonchiffon':'#FFFACD',
'lightblue':'#ADD8E6',
'lightcoral':'#F08080',
'lightcyan':'#E0FFFF',
'lightgoldenrodyellow':'#FAFAD2',
'lightgreen':'#90EE90',
'lightgray':'#D3D3D3',
'lightpink':'#FFB6C1',
'lightsalmon':'#FFA07A',
'lightseagreen':'#20B2AA',
'lightskyblue':'#87CEFA',
'lightslategray':'#778899',
'lightsteelblue':'#B0C4DE',
'lightyellow':'#FFFFE0',
'lime':'#00FF00',
'limegreen':'#32CD32',
'linen':'#FAF0E6',
'magenta':'#FF00FF',
'maroon':'#800000',
'mediumaquamarine':'#66CDAA',
'mediumblue':'#0000CD',
'mediumorchid':'#BA55D3',
'mediumpurple':'#9370DB',
'mediumseagreen':'#3CB371',
'mediumslateblue':'#7B68EE',
'mediumspringgreen':'#00FA9A',
'mediumturquoise':'#48D1CC',
'mediumvioletred':'#C71585',
'midnightblue':'#191970',
'mintcream':'#F5FFFA',
'mistyrose':'#FFE4E1',
'moccasin':'#FFE4B5',
'navajowhite':'#FFDEAD',
'navy':'#000080',
'oldlace':'#FDF5E6',
'olive':'#808000',
'olivedrab':'#6B8E23',
'orange':'#FFA500',
'orangered':'#FF4500',
'orchid':'#DA70D6',
'palegoldenrod':'#EEE8AA',
'palegreen':'#98FB98',
'paleturquoise':'#AFEEEE',
'palevioletred':'#DB7093',
'papayawhip':'#FFEFD5',
'peachpuff':'#FFDAB9',
'peru':'#CD853F',
'pink':'#FFC0CB',
'plum':'#DDA0DD',
'powderblue':'#B0E0E6',
'purple':'#800080',
'red':'#FF0000',
'rosybrown':'#BC8F8F',
'royalblue':'#4169E1',
'saddlebrown':'#8B4513',
'salmon':'#FA8072',
'sandybrown':'#FAA460',
'seagreen':'#2E8B57',
'seashell':'#FFF5EE',
'sienna':'#A0522D',
'silver':'#C0C0C0',
'skyblue':'#87CEEB',
'slateblue':'#6A5ACD',
'slategray':'#708090',
'snow':'#FFFAFA',
'springgreen':'#00FF7F',
'steelblue':'#4682B4',
'tan':'#D2B48C',
'teal':'#008080',
'thistle':'#D8BFD8',
'tomato':'#FF6347',
'turquoise':'#40E0D0',
'violet':'#EE82EE',
'wheat':'#F5DEB3',
'white':'#FFFFFF',
'whitesmoke':'#F5F5F5',
'yellow':'#FFFF00',
'yellowgreen':'#9ACD32'}

def matproj_formation_energy(formula):
    prop=['unit_cell_formula','material_id','energy_per_atom','formation_energy_per_atom']
    get_data=[]
    flag=False
    with open("mpb","r") as f:
        for i in f.readlines():
            if flag:
                if "search " in i:
                    break
                else:
                    get_data.append(json.loads(i.replace("\n","")))
            elif i.replace("search ","").replace("\n","")==formula:
                flag=True
    
    if not get_data:
        with open("mpb","a") as f:
            f.write("search "+formula+"\n")
            with MPRester("VUswQqBEWe4VFBZD25") as m:
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
            num=num*10+float(i)
        else:
            print(i)
    
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

def readtri(tri,formula_tuple,formation_energy):
    return [i*j[1]/formation_energy if not type(i)==bool else i for i,j in zip(tri,formula_tuple)]

def printtri(tri,formula_tuple,formation_energy):
    return [round(i/j[1]*formation_energy,3) for i,j in zip(tri,formula_tuple)]

def line2point(point1,point2):
    if abs(point1[0]-point2[0])<1e-2:
        temp=[(point1[1]-point2[1])/1e-2,(point1[0]*point2[1]-point2[0]*point1[1])/1e-2]
    else:
        temp=[(point1[1]-point2[1])/(point1[0]-point2[0]),(point1[0]*point2[1]-point2[0]*point1[1])/(point1[0]-point2[0])]
    return temp

def plot3(formula):
    #"Cs2SnCl6"
    formula_tuple=formula2tuple(formula)
    #(('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0))
    
    cp3=matproj_formation_energy("-".join([i[0] for i in formula_tuple]))
    #{(('Cl', 3.0), ('Cs', 1.0), ('Sn', 1.0)): ['mp-27394', -3.533219568, -1.8115027265862065], 
    #(('Cl', 5.0), ('Cs', 1.0), ('Sn', 2.0)): ['mp-30164', -3.551530498125, -1.6656268253663797], 
    #(('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0)): ['mp-608555', -3.3898603444444446, -1.922349615651341]}
    formation_energy=sum(list(zip(*formula_tuple))[1])*cp3[formula_tuple][-1]
    
    cp2={}
    
    for i in range(len(formula_tuple)):
        for j in range(len(formula_tuple))[i+1:]:
            cp2.update(matproj_formation_energy(formula_tuple[i][0]+"-"+formula_tuple[j][0]))
    # {'CsSn': ['mp-571056', -2.760116853125, -0.3076297395905172],
    #  'Cs4Sn23': ['mp-2496', -3.7147001966666666, -0.16648678358237556],
    #  'CsSn3': ['mp-865565', -3.38292671375, -0.15189722948275897],
    #  'CsCl': ['mp-573697', -3.31991964, -2.2549501239655174],
    #  'SnCl': ['mp-978928', -3.4328984175, -0.8108441600000003],
    #  'SnCl2': ['mp-569152', -3.561809344166667, -1.4022609525],
    #  'SnCl4': ['mp-29866', -3.2167423900000003, -1.4271986910000003]}
    
    for i,l in zip(cp3,cnames):
        if i==formula_tuple:
            # for j,k in zip(i,abc):
            #     print(k,j[0],printtri(xytotri(k),i,formation_energy))
            #     plt.annotate(j[0],k,xytext=(0,0),textcoords='offset points',fontsize=14)
            
            lines=[]
            for j in range(len(formula_tuple)):
                tri=[formation_energy/i[k][1] if j==k else 0 for k in range(len(formula_tuple))]
                xy=tritoxy(readtri(tri,formula_tuple,formation_energy))
                plt.annotate(formula_tuple[j][0],xy,xytext=(0,0),textcoords='offset points',fontsize=14)
                lines.append(xy)
            
            for j in range(len(formula_tuple)):
                for k in range(len(formula_tuple))[j+1:]:
                    temp=list(zip(lines[j],lines[k]))
                    plt.plot(temp[0],temp[1],color="k")
        
        else:
            temp_formula_num=np.array(list(zip(*i))[1])
            temp_formation_energy=sum(temp_formula_num)*cp3[i][-1]
            print(i,temp_formation_energy,l)
            for j in np.linspace(0,1,11):
                for k in np.linspace(0,1,11):
                    if not j+k>1:
                        xy=tritoxy([j,k,False])
                        pri=printtri(xytotri(xy),formula_tuple,formation_energy)
                        if np.array(pri).dot(temp_formula_num)>temp_formation_energy:
                            if 0 in pri:
                                print(pri)
                            plt.scatter(*xy,color=l)

    lines={}
    for i in cp2:
        tri=[]
        temp_formula=list(zip(*i))
        temp_formula_element=temp_formula[0]
        #temp_formula_num=np.array(temp_formula[1])
        temp_formation_energy=sum(temp_formula[1])*cp2[i][-1]
        for j in i:
            tri.append([])
            for k in formula_tuple:
                if k[0] in temp_formula_element:
                    if j[0]==k[0]:
                        print(j)
                        tri[-1].append(temp_formation_energy/j[1])
                    else:  
                        tri[-1].append(0)
                else:
                    tri[-1].append(False)
        
        temp0=tritoxy(readtri(tri[0],formula_tuple,formation_energy))
        temp1=tritoxy(readtri(tri[1],formula_tuple,formation_energy))
        lines[i]=line2point(temp0,temp1)
        plt.plot([temp0[0],temp1[0]],[temp0[1],temp1[1]])
        plt.annotate(tuple2formula(i),xy=temp0,xytext=(0,0),textcoords='offset points',fontsize=4)
        plt.annotate(tuple2formula(i),xy=temp1,xytext=(0,0),textcoords='offset points',fontsize=4)
    
    for i,j in zip(formula_tuple,labc):
        lines[i[0]+"_rich"]=j
    
    lines_set=[]
    for i in lines:
        for j in lines:
            if not i==j and not {i,j} in lines_set:
                print(i,j,printtri(xytotri(crosspoint(lines[i],lines[j])),formula_tuple,formation_energy))
                lines_set.append({i,j})
    
    ax=plt.gca()
    ax.set_aspect(1)
    #plt.savefig("pot.png",dpi=300)
    plt.show()
    return None

plot3("Cs2SnCl6")
