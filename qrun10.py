# -*- coding: utf-8 -*-
"""
Created on Mon May 31 21:38:56 2021

@author: dugue

version=0.1
"""

# !/home/phys/qif/anaconda3/bin/python

charge_corr_dict={'0Cs2AgInCl6': -0.10740636477538464,
'0Cs2HfCl6': -0.22865052200616062,
'0Cs2KInCl6c2h': -0.06571622202097328,
'0Cs2KInCl6oh': -0.15055106860460557,
#'0Cs2NaBiCl6': -0.1539912592779122,
'0Cs2NaBiCl6': -0.09500446458629402,
'0Cs2NaInCl6': -0.1968077259230447,
'0Cs2NaLaCl6': -0.15737025352611983,
'0Cs2NaScCl6': -0.1910407362901203,
'0Cs2NaYCl6': -0.1772608238927332,
'0Cs2SnCl6': -0.2562074178635485,
'0Cs2ZrCl6': -0.22401592892112798,
'0KCl': -0.36247277151921625/2,
'0NaCl': -0.4186868438803664/2,
'0Rb3InCl6c2h': -0.08056157539366501,
'0Rb3InCl6oh': -0.13942761574705112}

import numpy as np
import matplotlib.pyplot as plt

# python调用并回显linux命令
def linux_command(command):
    import os
    print(command)
    with os.popen(command, "r") as p:
        command_return=p.readlines()
    return command_return

import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
import numpy as np
import json
import random
import os

def matproj_formation_energy(formula):
    prop=['unit_cell_formula','material_id','energy_per_atom','formation_energy_per_atom']
    get_data=False
    with open("~/mpb","a+") as f:
        f.seek(0)
        for i in f.readlines():
            if type(get_data)==list:
                if "search " in i:
                    break
                else:
                    get_data.append(json.loads(i.replace("\n","")))
            elif i.replace("search ","").replace("\n","")==formula:
                get_data=[]
    
    if type(get_data)==bool:
        get_data=[]
        with open("mpb","a") as f:
            with MPRester("VUswQqBEWe4VFBZD25") as m:
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

abc=[[1/3**0.5,1],[0,0],[2/3**0.5,0]]
labc=[line2point(*i) for i in [[j for j in abc if not j==i] for i in abc]]

def parallel_line(ref_line,distance):
    if ref_line==labc[0]:
        return [ref_line[0],ref_line[1]+distance*(ref_line[0]**2+1)**0.5]
    else:
        return [ref_line[0],ref_line[1]-distance*(ref_line[0]**2+1)**0.5]

def distancetoline(line,xy):
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

# SOC crystal and mocular
soc_corr={"Bi":[-0.473,-0.808],
      "Cs":[-0.131,-0.143],
      "Na":[-0.000,-0.000],
      "Cl":[-0.007,-0.002]}

energy_corr={"Cl":0.614,"Bi":-0.808,"O":0.687}


def chempot_sample(formula,fixed={},start_point=[],other_element=[],soc=False,sample_times=100):
    with open(formula+"_".join([i+"_"+str(round(fixed[i],3)) for i in fixed])+".csv","w") as f:
        # "Cs2SnCl6" (('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0))
        formula_tuple=formula2tuple(formula)
        formation_mat={}
        other_formation_mat={}
        energy_per_atom={}
        for i in range(2**(len(formula_tuple)+len(other_element)))[1:]:
            bicode=bin(i)
            bicode=bicode.replace("0b","".join(["0" for j in range(len(formula_tuple)+len(other_element)-len(bicode)+2)]))
            ions=[]
            for j,k in zip(bicode,list(dict(formula_tuple).keys())+other_element):
                if j=="1":
                    ions.append(k)
            formation_dict=matproj_formation_energy("-".join(ions))
            # {'unit_cell_formula':['material_id','energy_per_atom','formation_energy_per_atom']
            # {(('Cl', 3.0), ('Cs', 1.0), ('Sn', 1.0)): ['mp-27394', -3.533219568, -1.8115027265862065], 
            # (('Cl', 5.0), ('Cs', 1.0), ('Sn', 2.0)): ['mp-30164', -3.551530498125, -1.6656268253663797],
            # (('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0)): ['mp-608555', -3.3898603444444446, -1.922349615651341]}
            
            for j in formation_dict:
                if len(j)==1:
                    energy_per_atom[j[0][0]]=round(formation_dict[j][1]+(energy_corr[j[0][0]] if j[0][0] in energy_corr else 0),3)
                else:
                    temp_formation_mat=[]
                    formation_energy=sum(list(zip(*j))[1])*formation_dict[j][-1]
                    if soc:
                        for k in dict(j):
                            if k in soc_corr:
                                formation_energy+=dict(j)[k]*(soc_corr[k][0]-soc_corr[k][1])
                    
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
                        other_formation_mat[j]=temp_formation_mat
                    else:
                        formation_mat[j]=temp_formation_mat
        
        print("element energy(soc&anion corr)")
        for i in energy_per_atom:
            print(f"{i:7} {energy_per_atom[i]:5.3f}")
        
        print("material",end="       ")
        real_formula_tuple=tuple([i for i in formula_tuple if not i[0] in fixed])
        for i in real_formula_tuple:
            print(i[0],end="  ")
        print("formation")
        
        for i in formation_mat:
            temp=f"{formation_mat[i][-1]:12}"
            for j in formation_mat[i][:-2]:
                temp+=f"{j:4.0f}"
            temp+=f"{formation_mat[i][-2]:8.2f}"
            print(temp)
        
        if other_formation_mat:
            print("material",end="       ")
            for i in real_formula_tuple:
                print(i[0],end="  ")
            for i in other_element:
                print(i,end="  ")
            print("formation")
            
            for i in other_formation_mat:
                temp=f"{other_formation_mat[i][-1]:12}"
                for j in other_formation_mat[i][:-2]:
                    temp+=f"{j:4.0f}"
                temp+=f"{other_formation_mat[i][-2]:8.2f}"
                print(temp)
                    

        chempot_min=[formation_mat[formula_tuple][-2]/i for i in formation_mat[formula_tuple][:-2]]
        
        f.write(",".join([str(i[0]) for i in real_formula_tuple])+"\n")
        f.write(",".join([str(i) for i in chempot_min])+"\n")
            
        for i in formation_mat:
            f.write(",".join([str(j) for j in formation_mat[i]])+"\n")
        
        if len(real_formula_tuple)==3:
            line_mat={}
            ax1=plt.subplot(111,aspect=1)
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
                ax1.plot(*list(zip(*end_points)),label=line_mat[i][4])   

            for i,j in zip(formula_tuple,labc):
                lines[i[0]+"_rich"]=j
            
            lines_set=[]
            for i in lines:
                for j in lines:
                    if not i==j and not {i,j} in lines_set:
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
            if len(real_formula_tuple)==3:
                ax1.scatter(*tritoxy(readtri(chempot_list,real_formula_tuple,formation_mat[formula_tuple][-2])),color="red")
            f.write("sample\n"+",".join([str(i) for i in chempot_list])+"\n")
            chempot_range={}
            for i in real_formula_tuple:
                chempot_range[(i[0],"max")]=chempot_list
                chempot_range[(i[0],"min")]=chempot_list
            
            ratio_list=[i/j for i,j in zip(chempot_list,chempot_min)]
            count=0
            while count<sample_times:
                temp_ratio_list=[i+0.1*(2*random.random()-1) for i in ratio_list[:-1]]
                temp_ratio_list.append(1-sum(temp_ratio_list))
                for i in temp_ratio_list:
                    if i>1 or i<0:
                        temp_ratio_list=[]
                        break
                    else:
                        chempot_list=[j*k for j,k in zip(temp_ratio_list,chempot_min)]
                
                if not temp_ratio_list:
                    continue
                
                for i in formation_mat:
                    if not i==formula_tuple:
                        if not condition(chempot_list,formation_mat[i][:-1]):
                            chempot_list=[]
                            break
                    
                if chempot_list:
                    count+=1
                    ratio_list=temp_ratio_list
                    f.write(",".join([str(i) for i in chempot_list])+"\n")
                    
                    for i in range(len(real_formula_tuple)):
                        chempot_range[(real_formula_tuple[i][0],"max")]=(chempot_list 
                                                                         if chempot_list[i]>chempot_range[(real_formula_tuple[i][0],"max")][i] 
                                                                         else chempot_range[(real_formula_tuple[i][0],"max")])
                        chempot_range[(real_formula_tuple[i][0],"min")]=(chempot_list 
                                                                         if chempot_list[i]<chempot_range[(real_formula_tuple[i][0],"max")][i] 
                                                                         else chempot_range[(real_formula_tuple[i][0],"min")])
            
            for i in chempot_range:
                chempot_range[i]=[round(j,3) for j in chempot_range[i]]
                
            result_dict={}
            for i in chempot_range:
                if tuple(chempot_range[i]) in result_dict:
                    result_dict[tuple(chempot_range[i])]["tag"].append("_".join(i))
                else:
                    result_dict[tuple(chempot_range[i])]={}
                    result_dict[tuple(chempot_range[i])]["tag"]=["_".join(i)]
                    result_dict[tuple(chempot_range[i])]["chempot"]={j[0]:k+energy_per_atom[j[0]] for j,k in zip(real_formula_tuple,chempot_range[i])}
            if other_element:
                for i in result_dict:
                    other_min=[]
                    for j in other_formation_mat.values():
                        other_min.append([j[-1],(j[-2]-sum([i[k]*j[k] for k in range(len(i))]))/j[-3]])
                    other_min.sort(key=lambda x:x[-1])
                    result_dict[i]["tag"].append(other_element[0]+"_"+other_min[0][0])
                    result_dict[i]["chempot"][other_element[0]]=round(other_min[0][1],3)+energy_per_atom[other_element[0]]
            if fixed:
                for i in result_dict:
                    for j in fixed:
                        result_dict[i]["chempot"][j]=round(fixed[j]+energy_per_atom[j],3)
            for i in result_dict:
                print(" ".join(result_dict[i]["tag"]))
                for j in result_dict[i]["chempot"]:
                    print(f"{j:5}{result_dict[i]['chempot'][j]:.3f}")
            if len(real_formula_tuple)==3:
                lines=[]
                for j,k in zip(real_formula_tuple,abc):
                    ax1.annotate(j[0]+"("+",".join([str(l) for l in printtri(xytotri(k),real_formula_tuple,formation_mat[formula_tuple][-2])])+")",
                                          k,xytext=(0,0),textcoords='offset points',fontsize=8)
                    lines.append(k)
               
                for j in range(3):
                    for k in range(3)[j+1:]:
                        temp=list(zip(lines[j],lines[k]))
                        ax1.plot(temp[0],temp[1],color="k")
                
                ax1.legend(fontsize=6)
                plt.savefig(formula+"_".join([i+"_"+str(round(fixed[i],3)) for i in fixed])+".png",dpi=300)
                plt.close('all')
        
        else:
            print("sample failed")
    return result_dict


# 不定分组匹配
def match_times(sep,char,string):
    # char="\w" for element "\d" for nelement
    import re
    # 保持格式一致
    if sep=="\s":
        string=" "+string
    else:
        string=sep+string
    match_result=re.search("("+sep+"+("+char+"+))",string)
    match_list=[]
    while match_result:
        unit=match_result.group(2)
        if char=="\d":
            unit=int(unit) 
        match_list.append(unit)
        match_result=re.search("("+sep+"+("+char+"+)){"+str(len(match_list)+1)+"}",string)
    return match_list

def csv_split(string):
    return string.replace("\n","").split(",,",1)[0].split(",")

def origin_name(distort_name):
    origin_name=False
    distort_name=distort_name.replace("_half","")
    if "p1" in distort_name:
        origin_name=distort_name.replace("p1","").split("_attempt")[0]
    elif "m1" in distort_name:
        origin_name=distort_name.replace("m1","").split("_attempt")[0]
    elif "ex" in distort_name:
        if "ex_ez" in distort_name:
            ex_name="ex_ez"
        elif "ex_cz" in distort_name:
            ex_name="ex_cz"
        elif "ex_cb" in distort_name:
            ex_name="ex_cb"
        elif "ex_vb" in distort_name:
            ex_name="ex_vb"
        elif "ex_nosym" in distort_name:
            ex_name="ex_nosym"
        else:
            ex_name="ex"
        origin_name=distort_name.replace(ex_name,"grd").split("_attempt")[0]
    return origin_name

def main_comp(spd_dict):
    max_value=0
    max_key=False
    for i in spd_dict:
        if spd_dict[i]>max_value:
            max_value=spd_dict[i]
            max_key=i
    return max_key
            
        
# 交互/非交互访问material project
def matproj(formula,prop,path):
    # prop_list=['band_gap','cif','density','diel','e_above_hull','elasticity',
    #            'elements','energy','energy_per_atom',
    #            'formation_energy_per_atom','full_formula','hubbards',
    #            'icsd_id','icsd_ids','is_compatible','is_hubbard',
    #            'material_id','nelements','nsites','oxide_type','piezo',
    #            'pretty_formula','spacegroup','tags','task_ids',
    #            'total_magnetization','unit_cell_formula','volume']
    prop_list=['band_gap','density','energy_per_atom',
               'formation_energy_per_atom','material_id','full_formula',
               'symbol','point_group','total_magnetization']
    point_group_dict={"1":"C1",
                      "4":"C4",
                      "2":"C2",
                      "2/m":"C2h",
                      "mm2":"C2v",
                      "3":"C3",
                      "-6":"C3h",
                      "-3":"C3i",
                      "3m":"C3v",
                      "4/m":"C4h",
                      "4mm":"C4v",
                      "6":"C6",
                      "6/m":"C6h",
                      "6mm":"C6v",
                      "-1":"Ci",
                      "m":"Cs",
                      "222":"D2",
                      "-42m":"D2d",
                      "mmm":"D2h",
                      "32":"D3",
                      "-3m":"D3d",
                      "-6m2":"D3h",
                      "422":"D4",
                      "4/mmm":"D4h",
                      "622":"D6",
                      "6/mmm":"D6h",
                      "432":"O",
                      "m-3m":"Oh",
                      "-4":"S4",
                      "23":"T",
                      "-43m":"Td",
                      "m-3":"Th"}
    from pymatgen.ext.matproj import MPRester
    with MPRester("VUswQqBEWe4VFBZD25") as m:
        if prop:
            m.get_structure_by_material_id(
                sorted(m.get_data(formula),key=lambda x:x[prop])[0]["material_id"],
                final=True,conventional_unit_cell=True).to("poscar",path+"/CONTCAR")
            return False
        else:
            for i in range(len(prop_list)):
                print(i,prop_list[i])
            key=input("proprty to show: ")
            if key=="all":
                key=prop_list
            else:
                key=[prop_list[i] for i in match_times("\s","\d",key)]
                
            value=[]
            for i in m.get_data(formula):
                value.append([])
                for j in key:
                    if j=='symbol':
                        value[-1].append(i['spacegroup'][j])
                    elif j=='point_group':
                        value[-1].append(point_group_dict[i['spacegroup'][j]])
                    else:
                        value[-1].append(i[j])
            value.sort(key=lambda x:x[key.index('formation_energy_per_atom')])
            
            key.insert(0,"# ")
            for i in key:
                print("%-10s" %i[:10],end=" ")
            print()
            for i in value:
                i.insert(0,value.index(i))
                for j in i:
                    if type(j)==float:
                        j=round(j,5)
                    print("%-10s" %j,end=" ")
                print()
            mp_list=match_times("\s","\d",input("structure to get: "))
            if len(mp_list)==1:
                m.get_structure_by_material_id(value[mp_list[0]][key.index("material_id")],
                                               final=True,conventional_unit_cell=True).to("poscar",path+"/CONTCAR")
                return False
            else:
                bundle_list=[]
                for i in mp_list:
                    linux_command("mkdir "+path+"/"+value[i][key.index("material_id")])
                    m.get_structure_by_material_id(value[i][key.index("material_id")],
                                                   final=True,conventional_unit_cell=True).to("poscar",
                                                                                              path+"/"+value[i][key.index("material_id")]+"/CONTCAR")
                    bundle_list.append(path+"/"+value[i][key.index("material_id")])
                return bundle_list

def plot_cc(energy_list,path):
    # energy_list 包含的能量依次为grd,ex,ex_at_grd,grd_at_ex
    import numpy as np
    import matplotlib.pyplot as plt
    
    c1=1
    c2=2
    num_size=12
    word_size=12
    
    grd,ex,ex_at_grd,grd_at_ex=energy_list
    c3=3
    plt.plot((c1,c1),(grd,ex_at_grd),'ro--')
    plt.plot((c1,c2),(ex_at_grd,ex),'ro--')
    plt.plot((c2,c2),(ex,grd_at_ex),'ro--')
    plt.plot((c2,c1),(grd_at_ex,grd),'ro--')
    plt.plot((c1,c2),(grd,ex),'ro--')
    
    plt.annotate("grd",xy=(c1,grd),xytext=(-35,-5),textcoords='offset points',fontsize=word_size)
    plt.annotate("ex@grd",xy=(c1,ex_at_grd),xytext=(-60,-5),textcoords='offset points',fontsize=word_size)
    plt.annotate("ex",xy=(c2,ex),xytext=(10,-5),textcoords='offset points',fontsize=word_size)
    plt.annotate("grd@ex",xy=(c2,grd_at_ex),xytext=(10,-5),textcoords='offset points',fontsize=word_size)
    
    # 标记电荷转变能
    plt.annotate(str(round(ex_at_grd-grd,2))+"eV",xy=(c1,(grd+ex_at_grd)/2),xytext=(-55,-5),textcoords='offset points',fontsize=num_size)
    plt.annotate(str(round(ex_at_grd-ex,2))+"eV",xy=((c1+c2)/2,(ex+ex_at_grd)/2),xytext=(10,20),textcoords='offset points',fontsize=num_size,
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    plt.annotate(str(round(ex-grd_at_ex,2))+"eV",xy=(c2,(ex+grd_at_ex)/2),xytext=(10,-5),textcoords='offset points',fontsize=num_size)
    plt.annotate(str(round(grd_at_ex-grd,2))+"eV",xy=((c1+c2)/2,(grd+grd_at_ex)/2),xytext=(10,-30),textcoords='offset points',fontsize=num_size,
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    plt.annotate(str(round(ex-grd,2))+"eV",xy=((c1+c2)/2,(grd+ex)/2),xytext=(-20,-5),textcoords='offset points',fontsize=num_size)
    
    x=np.linspace(2*c1-c2,3*c2-2*c1,100)
    coef=np.polyfit((c1,c2,c2+1e-5),(ex_at_grd,ex,ex),2)
    plt.plot(x,np.poly1d(coef)(x))
    
    x=np.linspace(3*c1-2*c2,2*c2-c1,10)
    coef=np.polyfit((c1,c1+1e-5,c2),(grd,grd,grd_at_ex),2)
    plt.plot(x,np.poly1d(coef)(x))
    
    plt.xlim(c1-3*(c2-c1),c2+3*(c2-c1))
    plt.xticks([])
    plt.xlabel('configurations')
    plt.ylim(grd-.6*(ex-grd),ex+.6*(ex-grd))
    plt.ylabel('energy(eV)')
    plt.title("configuration coordinate curves of "+path)
    # plt.legend()
    plt.savefig(path+".png",dpi=300)
    plt.close('all')
    return None

replace_dict={"ex_at_grd_pbe0":"sp@grd_w/o_soc",
"ex_at_grd_pbe0_soc":"sp@grd",
"ex_cb_at_grd_pbe0_soc":"cb@grd",
"ex_vb_at_grd_pbe0_soc":"vb@grd",
"ex_ez_pbe0_soc_temp1":"ez_static",
"grd_at_ex_cz_pbe0_soc":"grd@cz",
"grd_at_ex_ez_pbe0":"grd@ez_w/o_soc",
"grd_at_ex_ez_pbe0_soc":"grd@ez",
"grd_at_ex_cb_pbe0_soc":"grd@cb",
"grd_at_ex_vb_pbe0_soc":"grd@vb",
"relax_ex_cb_pbe0_soc":"cb",
"relax_ex_cz_pbe0_soc":"cz",
"relax_ex_ez_pbe0":"ez_w/o_soc",
"relax_ex_ez_pbe0_soc":"ez",
"relax_ex_vb_pbe0_soc":"vb",
"relax_grd_pbe0":"grd_w/o_soc",
"relax_grd_pbe0_soc":"grd",
"relax_grdm1_pbe0_soc":"1+_+_cbm",
"relax_grdp1_pbe0_soc":"1-_-_vbm",
"relax_ex_cz_half_pbe0_soc":"cz_half",
"relax_ex_ez_half_pbe0_soc":"ez_half",
"ex_vb_at_grd_m1_pbe0_soc":"1vb@grd",
"ex_vc_at_grd_m1_pbe0_soc":"1vc@grd",
"ex_at_grd_m1_pbe0_soc":"1sp@grd",
"ex_cb_at_grd_m1_pbe0_soc":"1cb@grd"}

def plot_cc2(energy_dict,path):
    # energy_list 包含的能量依次为grd,ex,ex_at_grd,grd_at_ex
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    word_size=14
    
    # line_dict={"relax_grd":(("relax","ex_at"),("relax","ex_cb_at"),("relax","ex_vb_at"),("relax_grd","ex_vb_at_grd_m1"),("relax_grd","ex_cb_at_grd_m1"),("relax_grd","ex_at_grd_m1"),
    #                         ("grd","ex_ez"),("grd","ex_cz"),("grd","ex_cb"),("grd","ex_vb"),
    #                         ("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #            "ex_at_grd":(("ex_at_grd","relax_ex_ez"),("ex_at_grd","relax_ex_cz")),
    #            "ex_cb_at_grd":(("ex_cb_at_grd","relax_ex_cb"),()),
    #            "ex_vb_at_grd":(("ex_vb_at_grd","relax_ex_vb"),()),
    #            "relax_ex":(("relax","grd_at"),()),
    #            "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1"),()),
    #            "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
    #            "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    
    # fit_dict={"relax_grd":(("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #           "relax_ex":(("relax_ex_ez","ex_at_grd"),("relax_ex_cz","ex_at_grd")),
    #           "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd"),("relax_ex_cb","ex_cb_at_grd_m1")),
    #           "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd"),("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
    #           "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    for i in ["relax_grd_pbe0","relax_ex_ez_pbe0","grd_at_ex_ez_pbe0","ex_at_grd_pbe0","ex_ez_pbe0_soc_temp1","ex_vb_at_grd_m1_pbe0_soc","ex_vc_at_grd_m1_pbe0_soc","ex_at_grd_m1_pbe0_soc","ex_cb_at_grd_m1_pbe0_soc"]:
        if i in energy_dict:
            del energy_dict[i]
    line_dict={"relax_grd":(("relax","ex_at","green"),("relax","ex_cb_at","green"),("relax","ex_vb_at","green"),("relax_grd","ex_vb_at_grd_m1"),("relax_grd","ex_cb_at_grd_m1"),("relax_grd","ex_at_grd_m1"),
                            ("grd","ex_ez"),("grd","ex_cz"),("grd","ex_cb"),("grd","ex_vb"),
                            ("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
                "ex_at_grd":(("ex_at_grd","relax_ex_ez"),("ex_at_grd","relax_ex_cz")),
                "ex_cb_at_grd":(("ex_cb_at_grd","relax_ex_cb"),()),
                "ex_vb_at_grd":(("ex_vb_at_grd","relax_ex_vb"),()),
                # "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1","red"),()),
                # "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
                # "relax_ex":(("relax","grd_at","blue"),()),
                "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1"),("relax","grd_at","red"),()),
                "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1"),("relax","grd_at","purple")),
                "relax_ex_ez":(("relax","grd_at","blue"),()),
                "relax_ex_cz":(("relax","grd_at","cyan"),()),
                "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    
    
    fit_dict={"relax_grd":(("relax_grd","grd_at_ex_ez","green"),("relax_grd","grd_at_ex_cz","green"),("relax_grd","grd_at_ex_cb","green"),("relax_grd","grd_at_ex_vb","green")),
              "relax_ex":(("relax_ex_ez","ex_at_grd","blue"),("relax_ex_cz","ex_at_grd","cyan")),
              "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd","red"),("relax_ex_cb","ex_cb_at_grd_m1")),
              "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd","purple"),("relax_ex_vb","ex_vb_at_grd_m1")),
              "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    ref=energy_dict["relax_grd_pbe0_soc"][1]
    xmax=-100
    xmin=100
    ymax=-100
    ymin=100
    for i in energy_dict:
        if type(energy_dict[i])==tuple:
            energy_dict[i]=(energy_dict[i][0],energy_dict[i][1]-ref)
            xmax=energy_dict[i][0] if energy_dict[i][0]>xmax else xmax
            xmin=energy_dict[i][0] if energy_dict[i][0]<xmin else xmin
            ymax=energy_dict[i][1] if energy_dict[i][1]>ymax else ymax
            ymin=energy_dict[i][1] if energy_dict[i][1]<ymin else ymin
    for i in energy_dict:
        if type(energy_dict[i])==tuple:
            for j in line_dict:
                if j in i:
                    for k in line_dict[j]:
                        if k:
                            if i.replace(*k[:2]) in energy_dict and type(energy_dict[i.replace(*k[:2])])==tuple:
                                pos=list(zip(energy_dict[i],energy_dict[i.replace(*k[:2])]))
                                if "relax_" in i and "relax_" in i.replace(*k[:2]):
                                    pass
                                else:
                                    if i.replace("relax_","")[:3]==i.replace(*k[:2]).replace("relax_","")[:3]:
                                        a=-10
                                    else:
                                        plt.plot(*pos,linestyle="dashed",marker="o",color=k[2] if len(k)==3 else "red")
                                        if "cb" in i or "cb" in i.replace(*k[:2]):
                                            a=-35
                                        else:
                                            a=5
                                        temp=str(round(abs(energy_dict[i.replace(*k[:2])][1]-energy_dict[i][1]),2))
                                        temp=temp if len(temp)==4 else temp+"0"
                                        plt.annotate(temp,
                                             xy=((pos[0][0]+pos[0][1])/2,(pos[1][0]+pos[1][1])/2),xytext=(a,10),
                                             textcoords='offset points',fontsize=12)
                                    # plt.annotate(str(round(abs(energy_dict[i.replace(*k[:2])][1]-energy_dict[i][1]),2)),
                                    #          xy=((pos[0][0]+pos[0][1])/2,(pos[1][0]+pos[1][1])/2),xytext=(a,10),
                                    #          textcoords='offset points',fontsize=12)
                            
            for j in fit_dict:
                if j in i:
                    for k in fit_dict[j]:
                        if k:
                            if i.replace(*k[:2]) in energy_dict and type(energy_dict[i.replace(*k[:2])])==tuple:
                                pos=list(zip(energy_dict[i],energy_dict[i.replace(*k[:2])]))
                                if len([l for l in energy_dict if "relax_ex_" in l])>1:
                                    if "relax_grd" in i:
                                    #if False:
                                        if pos[0][1]>pos[0][0]:
                                            x=np.linspace(pos[0][0]+2*abs(pos[0][0]-pos[0][1]),
                                                          pos[0][0]
                                                          ,50)
                                        else:
                                            x=np.linspace(pos[0][0],
                                                          pos[0][0]-1.2*abs(pos[0][0]-pos[0][1])
                                                          ,50)
                                    else:
                                        x=np.linspace(pos[0][0]-1.2*abs(pos[0][0]-pos[0][1]),
                                                      pos[0][0]+2*abs(pos[0][0]-pos[0][1]),100)
                                else:
                                    x=np.linspace(pos[0][0]-1.2*abs(pos[0][0]-pos[0][1]),
                                                  pos[0][0]+2*abs(pos[0][0]-pos[0][1]),100)
                                    
                                    
                                coef=np.polyfit((pos[0][0]-1e-5,pos[0][0]+1e-5,pos[0][1]),(pos[1][0],pos[1][0],pos[1][1]),2)
                                plt.plot(x,np.poly1d(coef)(x),color=k[2] if len(k)==3 else "red")
            #plt.plot(*energy_dict[i],"rx")
            if not "at" in i:
                if energy_dict["relax_grd_pbe0_soc"][0]>1:
                    if "cb" in i:
                        a=-30
                    elif "z" in i:
                        a=10
                    else:
                        a=-10
                else:
                    if "cz" in i:
                        a=-35
                    elif "z" in i:
                        a=0
                    else:
                        a=-10
                plt.annotate(replace_dict[i] if not replace_dict[i]=="grd" else "",xy=energy_dict[i],xytext=(a,30 if "relax_grd" in i else 10),textcoords='offset points',fontsize=14)
    plt.annotate(path.split("_")[-1],xy=(xmax,ymax+0.12*(ymax-ymin)),xytext=(0,0),textcoords='offset points',fontsize=14)
    from matplotlib.pyplot import MultipleLocator
    ax=plt.gca()
    if energy_dict["relax_grd_pbe0_soc"][0]>1:
        plt.xlabel('Average bond length ($\AA$)',fontsize=16)
        ax.xaxis.set_minor_locator(MultipleLocator(0.01))
    else:
        plt.xlabel('Distortion ($\AA$)',fontsize=16)
        ax.xaxis.set_minor_locator(MultipleLocator(0.02))
    plt.xlim(xmin-0.2*(xmax-xmin),xmax+1*(xmax-xmin))
    plt.ylim(ymin-0.1*(ymax-ymin),ymax+0.2*(ymax-ymin))
    plt.ylabel('Energy (eV)',fontsize=16)
    plt.tick_params(labelsize=14)
    #plt.title("configuration coordinate curves of "+path)
    # plt.legend()
    #plt.ylim(3.5,4.25)
    
    plt.tight_layout()
    
    
    ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    ax.set_aspect((xmax-xmin)/(ymax-ymin)*1.5)
    plt.savefig(path+".png",dpi=300)
    plt.close('all')
    return None


def plot_cc3(energy_dict,path,xmin,xmax,y1min,y1max,y2min,y2max,rat):
    # energy_list 包含的能量依次为grd,ex,ex_at_grd,grd_at_ex
    import numpy as np
    import matplotlib.pyplot as plt
    word_size=14
    
    # line_dict={"relax_grd":(("relax","ex_at"),("relax","ex_cb_at"),("relax","ex_vb_at"),("relax_grd","ex_vb_at_grd_m1"),("relax_grd","ex_cb_at_grd_m1"),("relax_grd","ex_at_grd_m1"),
    #                         ("grd","ex_ez"),("grd","ex_cz"),("grd","ex_cb"),("grd","ex_vb"),
    #                         ("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #            "ex_at_grd":(("ex_at_grd","relax_ex_ez"),("ex_at_grd","relax_ex_cz")),
    #            "ex_cb_at_grd":(("ex_cb_at_grd","relax_ex_cb"),()),
    #            "ex_vb_at_grd":(("ex_vb_at_grd","relax_ex_vb"),()),
    #            "relax_ex":(("relax","grd_at"),()),
    #            "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1"),()),
    #            "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
    #            "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    
    # fit_dict={"relax_grd":(("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #           "relax_ex":(("relax_ex_ez","ex_at_grd"),("relax_ex_cz","ex_at_grd")),
    #           "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd"),("relax_ex_cb","ex_cb_at_grd_m1")),
    #           "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd"),("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
    #           "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    for i in ["relax_grd_pbe0","relax_ex_ez_pbe0","grd_at_ex_ez_pbe0","ex_at_grd_pbe0","ex_ez_pbe0_soc_temp1","ex_vb_at_grd_m1_pbe0_soc","ex_vc_at_grd_m1_pbe0_soc","ex_at_grd_m1_pbe0_soc","ex_cb_at_grd_m1_pbe0_soc"]:
        if i in energy_dict:
            del energy_dict[i]
    print(energy_dict)
    # line_dict={"relax_grd":(("relax","ex_at"),("relax","ex_cb_at"),("relax","ex_vb_at"),("relax_grd","ex_vb_at_grd_m1"),("relax_grd","ex_cb_at_grd_m1"),("relax_grd","ex_at_grd_m1"),
    #                         ("grd","ex_ez"),("grd","ex_cz"),("grd","ex_cb"),("grd","ex_vb"),
    #                         ("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #            "ex_at_grd":(("ex_at_grd","relax_ex_ez"),("ex_at_grd","relax_ex_cz")),
    #            "ex_cb_at_grd":(("ex_cb_at_grd","relax_ex_cb"),()),
    #            "ex_vb_at_grd":(("ex_vb_at_grd","relax_ex_vb"),()),
    #            "relax_ex":(("relax","grd_at"),()),
    #            "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1"),()),
    #            "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
    #            "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    
    # fit_dict={"relax_grd":(("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #           "relax_ex":(("relax_ex_ez","ex_at_grd"),("relax_ex_cz","ex_at_grd")),
    #           "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd"),("relax_ex_cb","ex_cb_at_grd_m1")),
    #           "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd"),("relax_ex_vb","ex_vb_at_grd_m1")),
    #           "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    
    line_dict={"relax_grd":(("relax","ex_at","green"),("relax","ex_cb_at","green"),("relax","ex_vb_at","green"),("relax_grd","ex_vb_at_grd_m1"),("relax_grd","ex_cb_at_grd_m1"),("relax_grd","ex_at_grd_m1"),
                        ("grd","ex_ez"),("grd","ex_cz"),("grd","ex_cb"),("grd","ex_vb"),
                        ("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
            "ex_at_grd":(("ex_at_grd","relax_ex_ez"),("ex_at_grd","relax_ex_cz")),
            "ex_cb_at_grd":(("ex_cb_at_grd","relax_ex_cb"),()),
            "ex_vb_at_grd":(("ex_vb_at_grd","relax_ex_vb"),()),
            # "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1","red"),()),
            # "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
            # "relax_ex":(("relax","grd_at","blue"),()),
            "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1"),("relax","grd_at","red"),()),
            "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1"),("relax","grd_at","purple")),
            "relax_ex_ez":(("relax","grd_at","blue"),()),
            "relax_ex_cz":(("relax","grd_at","cyan"),()),
            "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    
    
    fit_dict={"relax_grd":(("relax_grd","grd_at_ex_ez","green"),("relax_grd","grd_at_ex_cz","green"),("relax_grd","grd_at_ex_cb","green"),("relax_grd","grd_at_ex_vb","green")),
              "relax_ex":(("relax_ex_ez","ex_at_grd","blue"),("relax_ex_cz","ex_at_grd","cyan")),
              "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd","red"),("relax_ex_cb","ex_cb_at_grd_m1")),
              "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd","purple"),("relax_ex_vb","ex_vb_at_grd_m1")),
              "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    ref=energy_dict["relax_grd_pbe0_soc"][1]
    
    f,(ax,ax2)=plt.subplots(2,1,sharex=True)
    
    for i in energy_dict:
        if type(energy_dict[i])==tuple:
            energy_dict[i]=(energy_dict[i][0],energy_dict[i][1]-ref)
    for i in energy_dict:
        if type(energy_dict[i])==tuple:
            for j in line_dict:
                if j in i:
                    for k in line_dict[j]:
                        if k:
                            if i.replace(*k[:2]) in energy_dict and type(energy_dict[i.replace(*k[:2])])==tuple:
                                pos=list(zip(energy_dict[i],energy_dict[i.replace(*k[:2])]))
                                if abs(pos[0][0]-pos[0][1])<1e-3:
                                    # or abs(pos[1][0]-pos[1][1])<3
                                    ax.plot(*pos,linestyle="dashed",marker="o",color=k[2] if len(k)==3 else "red")
                                    ax2.plot(*pos,linestyle="dashed",marker="o",color=k[2] if len(k)==3 else "red")
                                # else:
                                #     x=np.linspace(*pos[0],2)
                                #     print([x,-(x-pos[0][0])/(pos[0][0]-pos[0][1])*(pos[1][1]-pos[1][0]-2.92)])
                                #     ax2.plot(x,-(x-pos[0][0])/(pos[0][0]-pos[0][1])*(pos[1][1]-pos[1][0]-2.92),'ro--')
                                #     print([x,-(x-pos[0][1])/(pos[0][0]-pos[0][1])*(pos[1][1]-pos[1][0]-2.92)+pos[1][1]])
                                #     ax.plot(x,-(x-pos[0][1])/(pos[0][0]-pos[0][1])*(pos[1][1]-pos[1][0]-2.92)+pos[1][1],'ro--')
                                if y1max<(pos[1][0]+pos[1][1])/2<y2min:
                                    temp=str(round(abs(energy_dict[i.replace(*k[:2])][1]-energy_dict[i][1]),2))
                                    temp=temp if len(temp)==4 else temp+"0"
                                        
                                    if False:
                                        ax.annotate(temp,
                                                 xy=((pos[0][0]+pos[0][1])/2,((pos[1][0]+pos[1][1])/2-2)/8+y2min),xytext=(0,0),
                                                 textcoords='offset points',fontsize=word_size)
                                        ax2.annotate(temp,
                                                 xy=((pos[0][0]+pos[0][1])/2,((pos[1][0]+pos[1][1])/2-2)/8+y2min),xytext=(0,0),
                                                 textcoords='offset points',fontsize=word_size)
                                    else:
                                        ax.annotate(temp,
                                                 xy=((pos[0][0]+pos[0][1])/2,((pos[1][0]+pos[1][1])/2-2)/8+y1max),xytext=(0,0),
                                                 textcoords='offset points',fontsize=word_size)
                                        ax2.annotate(temp,
                                                 xy=((pos[0][0]+pos[0][1])/2,((pos[1][0]+pos[1][1])/2-2)/8+y1max),xytext=(0,0),
                                                 textcoords='offset points',fontsize=word_size)
                                        
                                # else:
                                #     ax.annotate(str(round(abs(energy_dict[i.replace(*k[:2])][1]-energy_dict[i][1]),2)),
                                #                  xy=((pos[0][0]+pos[0][1])/2,(pos[1][0]+pos[1][1])/2),xytext=(0,0),
                                #                  textcoords='offset points',fontsize=word_size)
                                #     ax2.annotate(str(round(abs(energy_dict[i.replace(*k[:2])][1]-energy_dict[i][1]),2)),
                                #                  xy=((pos[0][0]+pos[0][1])/2,(pos[1][0]+pos[1][1])/2),xytext=(0,0),
                                #                  textcoords='offset points',fontsize=word_size)
            for j in fit_dict:
                if j in i:
                    for k in fit_dict[j]:
                        if k:
                            if i.replace(*k[:2]) in energy_dict and type(energy_dict[i.replace(*k[:2])])==tuple:
                                pos=list(zip(energy_dict[i],energy_dict[i.replace(*k[:2])]))
                                if False:
                                    pass
                                # if "relax_grd" in i:
                                #     if pos[0][1]>pos[0][0]:
                                #         x=np.linspace(pos[0][0]+1.5*abs(pos[0][0]-pos[0][1]),
                                #                       pos[0][0]
                                #                       ,50)
                                #     else:
                                #         x=np.linspace(pos[0][0],
                                #                       pos[0][0]-1.5*abs(pos[0][0]-pos[0][1])
                                #                       ,50)
                                else:
                                    x=np.linspace(pos[0][0]-1.5*abs(pos[0][0]-pos[0][1]),
                                              pos[0][0]+1.5*abs(pos[0][0]-pos[0][1])
                                              ,100)
                                coef=np.polyfit((pos[0][0]-1e-5,pos[0][0]+1e-5,pos[0][1]),(pos[1][0],pos[1][0],pos[1][1]),2)
                                ax.plot(x,np.poly1d(coef)(x),color=k[2] if len(k)==3 else "red")
                                ax2.plot(x,np.poly1d(coef)(x),color=k[2] if len(k)==3 else "red")
            #ax.plot(*energy_dict[i],"rx")
            #ax2.plot(*energy_dict[i],"rx")
            if not "at" in i:
                ax.annotate("" if replace_dict[i]=="grd" else replace_dict[i],xy=energy_dict[i],xytext=(0,10),textcoords='offset points',fontsize=word_size)
                ax2.annotate("" if replace_dict[i]=="grd" else replace_dict[i],xy=energy_dict[i],xytext=(0,10),textcoords='offset points',fontsize=word_size)
    from matplotlib.pyplot import MultipleLocator
    #ax3=plt.gca()
    #ax.xaxis.set_major_locator(MultipleLocator(0.04))
    ax2.set_ylim(y1min,y1max)
    ax.set_ylim(y2min,y2max)
    ax.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop=False)
    ax2.xaxis.tick_bottom()
    d=.015
    plt.xlim(xmin,xmax)
    kwargs=dict(transform=ax.transAxes,color='k',clip_on=False)
    ax.plot((-d,+d),(-d,+d),**kwargs)
    ax.plot((1-d,1+d),(-d,+d),**kwargs)
    

    
    
    #ax2.annotate("Bi",xy=(2.78,4),xytext=(0,0),textcoords='offset points',fontsize=14)
    #ax.annotate("Bi",xy=(2.79,4),xytext=(0,0),textcoords='offset points',fontsize=14)
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d,+d),(1-d,1+d),**kwargs)
    ax2.plot((1-d,1+d),(1-d,1+d),**kwargs)
    #plt.subplots_adjust(left=0.125,bottom=0.1,right=0.9,top=0.9,wspace=0.2,hspace=0.35)
    plt.subplots_adjust(hspace=0.1)

    #plt.xlabel('Average bond length ($\AA$)',fontsize=16)
    plt.xlabel('Distortion ($\AA$)',fontsize=16)
    plt.ylabel('Energy (eV)',fontsize=16)
    ax.tick_params(labelsize=14)
    ax2.tick_params(labelsize=14)
    ax.set_aspect(rat)
    ax2.set_aspect(rat)
    plt.tight_layout()
    plt.savefig(path+".png",dpi=300)
    plt.close('all')
    return None

def plot_ctl(band_range,charge_list,path):
    # charge transition level
    import numpy as np
    import matplotlib.pyplot as plt
    
    x=np.linspace(band_range[0],band_range[1]+0.5,500)
    # 不能给出交点
    # plt.plot(x,[min([i[0]*j+i[1] for i in charge_list]) for j in x],color="red")
    
    y=[]
    for i in range(len(charge_list)):
        y.append(charge_list[i][0]*x+charge_list[i][1])
        
    node_list=[]
    y_min=[]
    prev_min_key=False
    for i in range(len(x)):
        min_key=False
        min_value=1e5
        for j in range(len(y)):
            if y[j][i]<min_value:
                min_value=y[j][i]
                min_key=j
        if type(prev_min_key)==bool:
            prev_min_key=min_key
        elif not min_key==prev_min_key:
            node_list.append((x[i],(y[min_key][i]+y[prev_min_key][i])*0.5))
            prev_min_key=min_key
        y_min.append(min_value)
        
    factor=2e-4
    
    y_min_max=max(y_min)
    y_min_min=min(y_min)
    
    plt.plot(x,y_min,color="red")
    plt.plot([band_range[0],band_range[1]],[y_min[0],y_min[-1]],"ro")
    plt.annotate(str(round(band_range[0],2))+"eV",xy=(band_range[0],y_min[0]),xytext=(-25,-15),textcoords='offset points',fontsize=14)
    plt.annotate(str(round(band_range[1],2))+"eV",xy=(band_range[1],y_min[-1]),xytext=(-20,15),textcoords='offset points',fontsize=14)
    plt.plot([band_range[0],band_range[0]],[y_min_min*(1+factor),y_min_max*(1-factor)],linestyle='dashed',color="blue")
    plt.plot([band_range[1],band_range[1]],[y_min_min*(1+factor),y_min_max*(1-factor)],linestyle='dashed',color="blue")
    
    for i in node_list:
        plt.plot(i[0],i[1],"ro",color="red")
        plt.annotate(str(round(i[0],2))+"eV",xy=(i[0],i[1]),xytext=(-20,10),textcoords='offset points',fontsize=14)
        plt.plot([i[0],i[0]],[y_min_min*(1+factor),y_min_max*(1-factor)],linestyle='dashed',color="blue")
    
    plt.xlabel('band(eV)')
    plt.ylim(y_min_min*(1+factor),)
    plt.ylabel('energy(eV)')
    plt.title("charge transition level of "+path)
    # plt.legend()
    plt.savefig(path+".png",dpi=300)
    plt.close('all')
    return list(zip(x,y))

def plot_ctl2(split_list,title,path):
    import matplotlib.pyplot as plt
    
    interval=1.5
    length=1
    

    # 价带顶对齐
    for i in split_list:
        avg=i[0]
        for j in range(len(i[:-1])):
            i[j]-=avg
    
    # 真空能级对齐
    # avg_dict={"KCl_Pb":1.614,
    #           "NaCl_Sn":1.999,
    #           "Cs2NaYCl6_Sb":2.312,
    #           "Cs2ZrCl6_Te":2.1470000000000002,
    #           "Cs2NaScCl6_As":2.38,
    #           "Cs2NaScCl6_Sb":2.38,
    #           "Cs2SnCl6_Te":3.2800000000000002,
    #           "Cs2ZrCl6_Se":2.1470000000000002,
    #           "Cs2NaLaCl6_Sb":2.508,
    #           "Cs2NaInCl6_Sb":2.497,
    #           "Cs2NaInCl6_Bi":2.497,
    #           "Cs2NaYBr6_Bi":2.184,
    #           "Cs2NaYCl6_Bi":2.312,
    #           "Cs2NaLaCl6_Bi":2.508}
    
    # for i in split_list:
    #     for j in range(len(i[:-1])):
    #         i[j]-=avg_dict[i[4]]
    
    line_list=[]
    for i in split_list:
        temp_line_list=[[i[0],[0],i[-1],[0,-40]],
                        [i[1],[1],str(round(i[1]-i[0],2))+"eV",[0,20]]]
        
        for j in i[2:-1]:
            temp_line_list.append([j,[i.index(j)],"",[0,-20]])
                
        line_list.append(temp_line_list)
        
    
    line_list.sort(key=lambda x:x[1][0]-x[0][0])
    line_list.sort(key=lambda x:x[0][2].split("/")[0].split("_")[-1])
    
    for i in range(len(line_list)):
        for j in range(len(line_list[i])):
            if j<2:
                plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="black")
            elif line_list[i][j][0]>500:
                line_list[i][j][0]-=1000
                plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],"r--")
            elif j==2:
                plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="red")
            else:
                plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="blue")
            if line_list[i][j][2]:
                plt.annotate(line_list[i][j][2],xy=(i*interval,line_list[i][j][0]),
                              xytext=line_list[i][j][3],textcoords='offset points',fontsize=6,rotation="270")
    plt.xticks([])
    plt.ylabel('energy(eV)')
    plt.title("ctl")
    if path:
        plt.savefig(path+".png",dpi=300)
    else:
        plt.show()
    plt.close("all")
    return None

def plot_ctl3(split_list,title,path):
    import matplotlib.pyplot as plt
    import random
    import numpy as np
    
    interval=1.5
    length=1
    

    # 价带顶对齐
    # for i in split_list:
    #     avg=i["VBM"]
    #     for j in i:
    #         if not j=="system":
    #             i[j]-=avg
    
    # 真空能级对齐
    avg_dict={"Cs2NaLaCl6":-5.584970968,
"Cs2SnCl6":-6.227959018,
"Cs2ZrCl6":-5.957924466,
"Cs2NaInCl6":-6.0492546,
"Cs2NaBiCl6":-5.844214162,
"KCl":-4.221406774,
"NaCl":-5.53,
"Cs2AgInCl6":-6.633276321,
"Cs2NaScCl6":-5.755407712,
"Cs2NaYCl6":-5.407844399,
"Rb3InCl6oh":-4.887988169,
"Rb3InCl6c2h":-4.887988169,
"Cs2HfCl6":-5.804700327,
"Cs2KInCl6oh":-5.515516321,
"Cs2KInCl6c2h":-5.515516321}
    
    for i in split_list:
        for j in i:
            if not j=="system":
                i[j]+=avg_dict[i["system"][1:].split("_")[0]]
    
    line_list=[]
    for i in split_list:
        temp_line_list=[]
        for j in i:
            if not j=="system":
                temp_line_list.append([i[j],[],j,[0,0]])
                # if j=="VBM":
                #     temp_line_list.append([i[j],[],i["system"],[0,-30]])
        line_list.append(temp_line_list)
        
    #line_list.sort(key=lambda x:max(list(zip(*x))[0]))
    for i in range(len(line_list)):
        for j in range(len(line_list[i])):
            if "VBM" in line_list[i][j][2]:
                plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="black",linewidth=1)
                plt.fill_between(np.linspace(i*interval-0.25,i*interval+length+0.25,15),line_list[i][j][0],-10,color="black",alpha=0.5)
            elif "CBM" in line_list[i][j][2]:
                plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="black",linewidth=1)
                plt.fill_between(np.linspace(i*interval-0.25,i*interval+length+0.25,15),line_list[i][j][0],2,color="black",alpha=0.25)
            elif line_list[i][j][2] in ["0/1","-1/0"]:
                if line_list[i][j][0]>1e4:
                    line_list[i][j][0]-=1e5
                    plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="red",linestyle="--",linewidth=1)
                else:
                    plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="red",linewidth=1)
            elif "Ag" in line_list[i][j][2]:
                plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="blue",linewidth=1)
            else:
                plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="cyan",linewidth=1)
                # if line_list[i][j][2]:
                #     if "_" in line_list[i][j][2] and not "/" in line_list[i][j][2]:
                #         plt.annotate(line_list[i][j][2],xy=(i*interval,line_list[i][j][0]),
                #                   xytext=line_list[i][j][3],textcoords='offset points',fontsize=4,rotation="330")
                #     else:
                #         plt.annotate(line_list[i][j][2],xy=(i*interval,line_list[i][j][0]),
                #                   xytext=line_list[i][j][3],textcoords='offset points',fontsize=4,rotation="0")
    plt.xticks(np.linspace(0,len(split_list),len(split_list)+1)[:-1]*interval+0.5,
               [i["system"][1:].replace("Cs2","").replace("Cl6","").replace("_Bi","").replace("_Sb","") for i in split_list])
    plt.ylabel('energy (eV)')
    plt.xlim(-0.25,14.75)
    plt.ylim(-9,1)
    #plt.title("ctl")
    if path:
        plt.savefig(path+".png",dpi=300)
    else:
        plt.show()
    plt.close("all")
    return None


    
def plot_ctl4(band_range,charge_list,path):
    # charge transition level
    # x=np.linspace(band_range[0]-0.2*(band_range[1]-band_range[0]),
    #               band_range[1]+0.5*(band_range[1]-band_range[0]),
    #               100)
    x=np.linspace(band_range[0],band_range[1],100)
    # 不能给出交点
    # plt.plot(x,[min([i[0]*j+i[1] for i in charge_list]) for j in x],color="red")
    
    for h in charge_list:
        y=[]
        for i in range(len(h)):
            y.append(h[i][0]*x+h[i][1])
            
        node_list=[]
        y_min=[]
        prev_min_key=False
        for i in range(len(x)):
            min_key=False
            min_value=1e5
            for j in range(len(y)):
                if y[j][i]<min_value:
                    min_value=y[j][i]
                    min_key=j
            if type(prev_min_key)==bool:
                prev_min_key=min_key
            elif not min_key==prev_min_key:
                node_list.append((x[i],(y[min_key][i]+y[prev_min_key][i])*0.5))
                prev_min_key=min_key
            y_min.append(min_value)
            
        factor=2e-4
        
        y_min_max=max(y_min)
        y_min_min=min(y_min)
        
        plt.plot(x,y_min,label=h[1][-1])
        plt.plot([band_range[0],band_range[1]],[y_min[0],y_min[-1]],"ro")
        plt.annotate(h[1][-1],xy=(band_range[0],y_min[0]),xytext=(0,0),textcoords='offset points',fontsize=6)
        #plt.annotate(h[1][-1],xy=(band_range[1],y_min[-1]),xytext=(-20,15),textcoords='offset points',fontsize=14)
        plt.annotate(h[1][-1],xy=(band_range[1],y_min[-1]),xytext=(0,0),textcoords='offset points',fontsize=6)
        plt.plot([band_range[0],band_range[0]],[y_min_min*(1+factor),y_min_max*(1-factor)],color="blue")
        plt.plot([band_range[1],band_range[1]],[y_min_min*(1+factor),y_min_max*(1-factor)],color="blue")
        
        for i in node_list:
            plt.plot(i[0],i[1],"ro")
            #plt.annotate(h[1][-1],xy=(i[0],i[1]),xytext=(0,0),textcoords='offset points',fontsize=6)
            #plt.plot([i[0],i[0]],[y_min_min*(1+factor),y_min_max*(1-factor)],linestyle='dashed',color="blue")
    
    plt.xlabel('band(eV)')
    plt.xlim(band_range[0]-0.1*(band_range[1]-band_range[0]),band_range[1]+0.3*(band_range[1]-band_range[0]))
    #plt.ylim(y_min_min*(1+factor),)
    plt.ylabel('energy(eV)')
    #plt.title("charge transition level of "+path)
    plt.legend(loc="upper right")
    #plt.show()
    plt.savefig(path+".png",dpi=300)
    plt.close('all')
    return list(zip(x,y))

def plot_ctl5(band_range, charge_list, temp,path):
    x = np.linspace(*band_range, 100)
    fermi=[]
    for i in charge_list:
        y = []
        for j in i:
            y.append([[j[0]*k+j[1],j[0]] for k in x])

        node_list = []
        y_min = []
        for j in range(len(x)):
            min_value = [1e5]
            for k in range(len(y)):
                if y[k][j][0] < min_value[0]:
                    min_value = y[k][j]
            if y_min and not min_value[1] == y_min[-1][1]:
                node_list.append((x[j], (min_value[0]+y_min[-1][0])*0.5))
            y_min.append(min_value)
        fermi.append(y_min)

        # y_min_max=list(zip(*y_min))[0]
        # y_min_min = min(y_min_max)
        # y_min_max = max(y_min_max)

        plt.plot(x, list(zip(*y_min))[0], label=i[1][-1],linestyle="dashed" if "o" in i[1][-1] else "solid")
        # plt.plot([band_range[0], band_range[1]], [y_min[0][0], y_min[-1][0]], "ro")
        #plt.annotate(i[1][-1], xy=(x[0], y_min[0][0]),xytext=(0, 0), textcoords='offset points', fontsize=6)
        #plt.annotate(i[1][-1], xy=(x[-1], y_min[-1][0]),xytext=(0, 0), textcoords='offset points', fontsize=6)
        # plt.plot([band_range[0], band_range[0]], [y_min_min *(1+factor), y_min_max*(1-factor)], color="blue")
        # plt.plot([band_range[1], band_range[1]], [y_min_min *(1+factor), y_min_max*(1-factor)], color="blue")

        # for j in node_list:
        #     plt.plot(j[0], j[1], "ro")
            #plt.annotate(i[1][-1], xy=(j[0], j[1]), xytext=(0, 0),textcoords='offset points', fontsize=6)
    
    fermi=list(zip(*fermi))
    y_min=1e5
    y_max=-1e5
    for i in fermi:
        atemp=list(zip(*i))[0]
        y_min=min(min(atemp),y_min)
        y_max=max(max(atemp),y_max)
    factor = 2e-4
    plt.plot([band_range[0], band_range[0]], [y_min *(1+factor), y_max*(1-factor)], color="blue")
    plt.plot([band_range[1], band_range[1]], [y_min *(1+factor), y_max*(1-factor)], color="blue")
    plt.annotate("  VBM", xy=(x[0], y_min),xytext=(0, 0), textcoords='offset points', fontsize=12)
    plt.annotate("  CBM", xy=(x[-1], y_min),xytext=(0, 0), textcoords='offset points', fontsize=12)
    for i,j in zip(fermi,x):
        #print(j,sum([k[1]*np.exp(-k[0]/temp) for k in i]))
        if sum([k[1]*np.exp(-k[0]/temp) for k in i])<0:
            plt.plot([j,j],[y_min,y_max], "b--")
            plt.annotate("  Fermi level="+str(round(j,2))+"eV under "+str(round(temp/0.025852*300))+"K", xy=(j, y_max), xytext=(-30, 0),textcoords='offset points', fontsize=12)
            break

    plt.xlabel('Band (eV)')
    plt.xlim(band_range[0]-0.1*(band_range[1]-band_range[0]),
             band_range[1]+0.3*(band_range[1]-band_range[0]))
    #plt.ylim(y_min_min*(1+factor),)
    plt.ylabel('Energy (eV)')
    plt.legend(loc="upper right")
    #plt.title("charge transition level of "+path)
    # plt.legend()
    #plt.show()
    plt.savefig(path+".png", dpi=300)
    plt.close('all')
    return list(zip(x, y))


def plot_split(split_list,path):
    import numpy as np
    import matplotlib.pyplot as plt
    
    interval=1.5
    length=1
    cfs,both,sos=split_list
    cfs=np.array(cfs)-sum(cfs)/6
    both=np.array(both)-sum(both)/6
    sos=np.array(sos)-sum(sos)/6
    line_list=[[[0,         [0,1,2,3,4,5],  "p orbitals",                                                                [-10,10]]],
                [[cfs[0],   [0],            str(round(cfs[1]-cfs[0],2))+"eV" if cfs[1]-cfs[0]>0.05 else "",     [0,50]],
                 [cfs[1],   [1],            str(round(cfs[2]-cfs[1],2))+"eV" if cfs[2]-cfs[1]>0.05 else "",     [0,50]],
                 [cfs[2],   [2],            str(round(cfs[3]-cfs[2],2))+"eV" if cfs[3]-cfs[2]>0.05 else "",     [0,50]],
                 [cfs[3],   [3],            str(round(cfs[4]-cfs[3],2))+"eV" if cfs[4]-cfs[3]>0.05 else "",     [0,50]],
                 [cfs[4],   [4],            str(round(cfs[5]-cfs[4],2))+"eV" if cfs[5]-cfs[4]>0.05 else "",     [0,50]],
                 [cfs[5],   [5],            "lf",                                                              [10,10]]],
                [[both[0],  [0],            str(round(both[1]-both[0],2))+"eV" if both[1]-both[0]>0.05 else "", [0,50]],
                 [both[1],  [1],            str(round(both[2]-both[1],2))+"eV" if both[2]-both[1]>0.05 else "", [0,50]],
                 [both[2],  [2],            str(round(both[3]-both[2],2))+"eV" if both[3]-both[2]>0.05 else "", [0,50]],
                 [both[3],  [3],            str(round(both[4]-both[3],2))+"eV" if both[4]-both[3]>0.05 else "", [0,50]],
                 [both[4],  [4],            str(round(both[5]-both[4],2))+"eV" if both[5]-both[4]>0.05 else "", [0,50]],
                 [both[5],  [5],            "lf + soc",                                                          [0,10]]],
                [[sos[0],   [0],            str(round(sos[1]-sos[0],2))+"eV" if sos[1]-sos[0]>0.05 else "",     [0,50]],
                 [sos[1],   [0],            str(round(sos[2]-sos[1],2))+"eV" if sos[2]-sos[1]>0.05 else "",     [0,50]],
                 [sos[2],   [0],            str(round(sos[3]-sos[2],2))+"eV" if sos[3]-sos[2]>0.05 else "",     [0,50]],
                 [sos[3],   [0],            str(round(sos[4]-sos[3],2))+"eV" if sos[4]-sos[3]>0.05 else "",     [0,50]],
                 [sos[4],   [0],            str(round(sos[5]-sos[4],2))+"eV" if sos[5]-sos[4]>0.05 else "",     [0,50]],
                 [sos[5],   [0],            "soc",                                                              [10,10]]],
                [[0,        [],             "",                                                                [10,10]]]]
    
    for i in range(len(line_list)):
        for j in range(len(line_list[i])):
            plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="red")
            if line_list[i][j][2] and not "eV" in line_list[i][j][2]:
                plt.annotate(line_list[i][j][2],xy=(i*interval,line_list[i][j][0]),xytext=line_list[i][j][3],textcoords='offset points',fontsize=14)
            if not i==len(line_list)-1:
                for k in line_list[i][j][1]:
                    plt.plot([i*interval+length,(i+1)*interval],[line_list[i][j][0],line_list[i+1][k][0]],color="red",linestyle='dashed')
    
    plt.xticks([])
    #plt.xlabel('                  ex-e          ex-e+soc        grd+soc        grd')
    plt.ylim(min([min(i) for i in [cfs,sos,both]])-0.25,max([max(i) for i in [cfs,sos,both]])+0.25)
    plt.ylabel('energy(eV)')
    plt.title("crystal field splitting and spin orbital splitting of "+path)
    plt.savefig(path+".png",dpi=300)
    plt.close("all")
    return None

def plot_split_cc(split_list,path):
    import matplotlib.pyplot as plt
    
    x=[i[1] for i in split_list]
    for i in range(len(split_list[0][0])):
        plt.plot(x,[j[0][i] for j in split_list],color="red")
    plt.ylabel('energy(eV)')
    plt.title("crystal field splitting cc"+path)
    plt.savefig(path+".png",dpi=300)
    plt.close("all")
    return None

def plot_split2(split_list,title,path):
    import matplotlib.pyplot as plt
    
    interval=1.5
    length=1
    
    for i in split_list:
        avg=sum(i[:6])/6
        for j in range(len(i[:6])):
            i[j]-=avg
    line_list=[]
    for i in split_list:
        line_list.append([[i[0],[0],i[6],[0,-40]],
                          [i[1],[1],"",[0,50]],
                          [i[2],[2],"",[0,50]],
                          [i[3],[3],"",[0,50]],
                          [i[4],[4],"",[0,50]],
                          [i[5],[5],str(round(i[5]-i[0],2))+"eV",[0,10]]])
    line_list.sort(key=lambda x:x[5][0]-x[0][0])
    
    for i in range(len(line_list)):
        for j in range(len(line_list[i])):
            plt.plot([i*interval,i*interval+length],[line_list[i][j][0],line_list[i][j][0]],color="red")
            if line_list[i][j][2]:
                plt.annotate(line_list[i][j][2],xy=(i*interval,line_list[i][j][0]),
                              xytext=line_list[i][j][3],textcoords='offset points',fontsize=6,rotation="270")
    
    plt.xticks([])
    plt.ylabel('energy(eV)')
    plt.title("ctl")
    if path:
        plt.savefig(path+".png",dpi=300)
    else:
        plt.show()
    plt.close("all")
    return None

def plot_nei(nll,path):
    import matplotlib.pyplot as plt
    
    fig,ax=plt.subplots(len(nll),sharex=True,sharey=True)
    for i in range(len(nll)):
        if i==0:
            ax[i].set_title("neighbours of "+path)
        ax[i].set_yticks([])
        ax[i].annotate(nll[i][-1],xy=(3,6),xytext=(80,-10),textcoords='offset points',fontsize=12)
        for j in nll[i][:-1]:
            ax[i].plot([j[1],j[1]],[0,j[2]],"r" if j[0]=="Cl" else "b",label=j[0])

    plt.subplots_adjust(hspace=0)
    plt.xlabel('distance(Å)')
    plt.savefig(path+".png",dpi=300)
    plt.close("all")
    return None

def plot_nei2(nll,title,path):
    # freq distribution
    # import numpy as np
    import matplotlib.pyplot as plt
    
    fig,ax=plt.subplots(len(nll),sharex=True,sharey=True)
    nll.sort(key=lambda x:x[0][1])
    for i in range(len(nll)):
        if i==0:
            ax[i].set_title("neighbours of "+path)
        ax[i].set_yticks([])
        ax[i].annotate(nll[i][-1],xy=(3,6),xytext=(80,-7),textcoords='offset points',fontsize=12)
        for j in nll[i][:-1]:
            ax[i].plot([j[1],j[1]],[0,j[2]],"r" if j[0]=="Cl" else "b",label=j[0])

    plt.subplots_adjust(hspace=0)
    plt.xlabel('distance(Å)')
    if path:
        plt.savefig(path+".png",dpi=300)
    else:
        plt.show()
    plt.close("all")
    return None
    
def plot_pot(x,y,path):
    import matplotlib.pyplot as plt
    plt.plot(x[1:],y[1:],'r')
    plt.xlabel(x[0])
    plt.ylabel(y[0])
    plt.title(path)
    plt.savefig(path+".png",dpi=300)
    plt.close('all')
    return None

def plot_band(band,label,path):
    import matplotlib.pyplot as plt
    for i in band[1:]:
        if max(i)>-10:
            plt.plot(band[0],i,"r")
    for i in label:
        plt.xticks(label[0],label[1])
        plt.plot([label[0],label[0]],[-10,10],"b--")
        
    # plt.xlabel(x[0])
    # plt.ylabel(y[0])
    # plt.title(path)
    plt.savefig(path+".png",dpi=300)
    plt.close('all')
    return None

def read_procar(path):
    import re
    band_pattern=re.compile("band\s+(\d+) # energy\s+([\d\.\-]+) # occ.\s+([\d\.]+)")
    p=poscar()
    p.read(path+"/CONTCAR")
    spd_mode=False
    spd_begin=False
    spd_dict={}
    current_band={}
    band_list=[]
    procar_list=[]
    
    with open(path+"/PROCAR","r") as f:
        for i in f.readlines():
            # band行
            if "occ." in i:
                current_band["num"],current_band["energy"],current_band["occ"]=[float(j) for j in band_pattern.search(i).groups()]
                if current_band["num"]==1 and band_list:
                    procar_list.append(band_list)
                    band_list=[]
            # ion行
            elif " s " in i:
                if not spd_mode:
                    spd_mode=match_times("\s","[\w\-]",i)
                    spd_pattern=re.compile("\s+".join(["([\d\.]+)" for j in range(len(spd_mode))]))
                spd_begin=True
            elif spd_begin:
                match_result=spd_pattern.search(i)
                if match_result:
                    # tot>0.01 获取主要部分
                    if float(match_result.groups()[-1])>0.01:
                        for j,k in zip(spd_mode[1:-1],match_result.groups()[1:-1]):
                            # 投影>0.001 获取主要部分
                            if float(k)>0.001:
                                component=(p.label[int(match_result.groups()[0])-1][0],j.replace("x","d")[0])
                                if not component in spd_dict:
                                    #spd_dict.append(component)
                                    spd_dict[component]=float(k)
                                else:
                                    spd_dict[component]+=float(k)
                else:
                    # 当前band结束
                    spd_begin=False
                    if spd_dict:
                        for j in spd_dict:
                            spd_dict[j]=round(spd_dict[j],3)
                    current_band["spd"]=spd_dict
                    spd_dict={}
                    band_list.append(current_band)
                    current_band={}
        procar_list.append(band_list)
    return procar_list

def energy(path):
    #import re
    with open(path+"/OSZICAR","r") as f:
        # match_result=re.search("E0=\s+([\+\-\.\dE]+)",f.readlines()[-1])
        # if match_result:
        #     return float(match_result.group(1))
        print(path,end=" ")
        last2lines=f.readlines()[-2:]
        if "E0" in last2lines[-1]:
            print(float(last2lines[-1].split("E0= ")[1].split(" ")[0]))
            return float(last2lines[-1].split("E0= ")[1].split(" ")[0])
        elif "E+" in last2lines[-1]:
            print(float([i for i in last2lines[-1].split(" ") if len(i)>1 and i[:2]=="-0"][0]))
            return float([i for i in last2lines[-1].split(" ") if len(i)>1 and i[:2]=="-0"][0])
        elif "E+" in last2lines[-2]:
            print(float([i for i in last2lines[-2].split(" ") if len(i)>1 and i[:2]=="-0"][0]))
            return float([i for i in last2lines[-2].split(" ") if len(i)>1 and i[:2]=="-0"][0])
        else:
            return False
        
class poscar:
    def __init__(self):
        pass
    
    def read(self,path):
        import numpy as np
        with open(path,"r") as f:
            lines=f.readlines()
        self.comment=lines[0]
        self.scale=float(lines[1].replace("\n",""))
        self.axis=np.array([[float(j) for j in lines[i].replace("\n","").split(" ") if j] for i in range(5)[2:]])
        self.element=[i for i in lines[5].replace("\n","").split(" ") if i]
        self.number=[int(i) for i in lines[6].replace("\n","").split(" ") if i]
        self.label=[]
        # in case 相同元素未在一起
        for i,j in zip(self.element,self.number):
            for k in range(j+1)[1:]:
                while (i,k) in self.label:
                    k+=1
                self.label.append((i,k))
        self.mode=lines[7]
        self.position=[]
        if "Direct configuration=" in lines[7]:
            temp_position=[]
            for i in range(len(lines))[8:]:
                if "Direct configuration=" in lines[i]:
                    self.position.append(temp_position)
                    temp_position=[]
                else:
                    temp_position.append(np.array([float(j) for j in lines[i].replace("\n","").split(" ")]))

        else:
            for i in range(len(lines))[8:]:
                temp_position=[]
                count=0
                for j in lines[i].replace("\n","").split(" "):
                    if j and count<3:
                        temp_position.append(float(j))
                        count+=1
                if temp_position:
                    self.position.append(np.array(temp_position))
                else:
                    break
            self.layer=[]
            layer=[]
            z=0
            for i in sorted([(i[0],i[1],j) for i,j in zip(self.label,self.position)],key=lambda x:x[2][2]):
                if (i[2][2]-z)<0.05:
                    layer.append(i)
                else:
                    self.layer.append(layer)
                    layer=[i]
                    z=i[2][2]
            self.layer.append(layer)
                
            
        return None
    
    
    def atom(self,element):
        import numpy as np
        self.comment=element+"\n"
        self.scale=1
        self.axis=np.array([[10,0,0],
                            [0,10,0],
                            [0,0,10]])
        self.element=[element]
        self.number=[1]
        self.mode="Direct\n"
        self.position=[np.array([0.5,0.5,0.5])]
        return None
    
    def write(self,path):
        with open(path,"w") as f:
            f.write(self.comment)
            f.write(str(self.scale)+"\n")
            for i in self.axis:
                for j in i:
                    f.write(str(j)+" ")
                f.write("\n")
            for j in self.element:
                f.write(str(j)+" ")
            f.write("\n")
            for j in self.number:
                f.write(str(j)+" ")
            f.write("\n")
            f.write(self.mode)
            for i in self.position:
                for j in i:
                    f.write(str(j)+" ")
                f.write("\n")
        return None
    
    def copy(self,full=False):
        p=poscar()
        p.comment=self.comment
        p.scale=self.scale
        p.axis=self.axis
        p.element=self.element[:]
        p.number=self.number[:]
        p.label=self.label[:]
        p.mode=self.mode
        if full:
            p.position=self.position[:]
        else:
            p.position=[]
        return p
        
    def __add__(self,other_poscar):
        result=self.copy()
        result.position=[i + j for i,j in zip(self.position, other_poscar.position)]
        return result
    
    def __sub__(self,other_poscar):
        result=self.copy()
        for i in range(len(self.position)):
            diff=self.position[i]-other_poscar.position[i]
            for j in range(len(diff)):
                if diff[j]>0.5:
                    diff[j]-=1
                elif diff[j]<-0.5:
                    diff[j]+=1
            result.position.append(diff)
        return result
    
    def __mul__(self,multiplier):
        result=self.copy()
        result.position=[i*multiplier for i in self.position]
        return result
    
    def move(self,vector):
        for i in range(len(self.position)):
            self.position[i]+=vector
            for j in range(3):
                if self.position[i][j]>=1:
                    self.position[i][j]-=1
                elif self.position[i][j]<0:
                    self.position[i][j]+=1
        return None
    
    # configuration coordinate curves
    def cc(self,vib,ini,fin,step):
        poscar_list=[]
        for i in range(step):
            factor=round(ini+i*(fin-ini)/(step-1),2)
            poscar_list.append((factor,self+vib*factor))
        return poscar_list
    
    def cc2(self):
        import numpy as np
        z2x=self.copy()
        z2y=self.copy()
        for i in self.position:
            z2x.position.append(np.array([i[2],i[0],i[1]]))
            z2y.position.append(np.array([i[1],i[2],i[0]]))
        temp_pos=[]
        for i,j in zip(z2x.position,z2x.label):
            temp_pos.append(z2y.position[z2y.label.index(z2y.distance(i,False)[0][0])])
        z2y.position=temp_pos
        return z2x-z2y
    
    def avg(self):
        self_avg={}
        for i,j in zip(self.label,self.position):
            temp_position=[]
            for k in range(3):
                if j[k]>0.9:
                    #temp_position.append(j[k]-0.5)
                    temp_position.append(j[k]-1)
                # elif j[k]<0.1:
                #     temp_position.append(j[k]+0.5)
                else:
                    temp_position.append(j[k])
                    
            temp_position=np.array(temp_position)
            print(temp_position)
            if i[0] in self_avg:
                self_avg[i[0]]=[self_avg[i[0]][0]+temp_position,self_avg[i[0]][1]+1]
            else:
                self_avg[i[0]]=[temp_position,1]
        
        return self_avg
    
    # def mimic(self,ref_poscar):
    
    def distort(self,mode):
        import numpy as np
        import random
        p=poscar()
        if mode=="random":
            p.position=[[np.array((random.random()-.5)*5e-3) for j in i]for i in self.position]
            return self+p
        elif mode=="ez":
            xyz=(-1,-1,2)
        elif mode=="cz":
            xyz=(1,1,-2)
        p.position=[np.array(
            [xyz[k]*1e-3*(np.sign(j[k]-0.5)-2*j[k]+1) if not 0.51>j[k]>0.49 else 0 for k in range(3)]
            ) for j in self.position]
        return self+p
    
    def distort_sphere(self,distort_list):
        import numpy as np
        p=self.copy(full=True)
        for i,j in zip(self.distance(self.label[0],False)[1:],distort_list):
            temp_position=[]
            for k in self.position[self.label.index(i[0])]-self.position[0]:
                if k>0.5:
                    k-=1
                elif k<-0.5:
                    k+=1
                k*=1+j/i[1]
                if k>0.5:
                    k-=1
                elif k<-0.5:
                    k+=1
                temp_position.append(k)
            p.position[self.label.index(i[0])]=np.array(temp_position)+self.position[0]
        return p
    
    def distance(self,site,full):
        import numpy as np
        if type(site)==tuple:
            site=self.position[self.label.index(site)]
        distance_list=[]
        for i,j in zip(self.label,self.position):
            avector=j-site
            # 距离
            if full:
                distance_list.append([i,np.round(np.array([i if abs(i)<0.5 else (i+1 if i<0 else i-1) for i in avector]).dot(self.axis),5)])
            else:
                distance_list.append([i,round(np.linalg.norm(np.array([i if i<0.5 else 1-i for i in abs(avector)]).dot(self.axis)),5)])
        distance_list.sort(key=lambda x:x[0])
        distance_list.sort(key=lambda x:np.linalg.norm(x[1]))
        return distance_list
    
    def cc_xtick(self,mode,num):
        distance_list=[]
        for i in self.distance(self.label[0],False)[1:]:
            if len(distance_list)<num:
                distance_list.append(i[1])
        
        avg=sum(distance_list)/num
        if mode=="avg":
            return avg
        else:
            if distance_list[0]-avg>-1e-3:
                return 0
            elif (distance_list[-1]-avg)/(distance_list[0]-avg)<-1:
                return distance_list[-1]-avg
            else:
                return distance_list[0]-avg
            
    def env(self,site,prec,full):
        distance_list=self.distance(site,False)
        for i in distance_list:
            i[1]=round(i[1],prec)
        env_list=[]
        if full:
            for i in distance_list[1:]:
                if env_list and i[0][0]==env_list[-1][0] and i[1]==env_list[-1][1]:
                    env_list[-1][2]+=1
                    env_list[-1][3].append(i[0])
                else:
                    env_list.append([i[0][0],i[1],1,[i[0]]])
        else:
            for i in distance_list[1:]:
                if env_list and i[0][0]==env_list[-1][0] and i[1]==env_list[-1][1]:
                    env_list[-1][2]+=1
                else:
                    env_list.append([i[0][0],i[1],1])
        return env_list
    
    def sites(self,site=False,nei=-1,prec=2):
        site_list=[]
        for i in self.label:
            if site and not site==i[0]:
                continue
            site_env=self.env(i,prec,False)[:nei]
            if site_list:
                site_saved=False
                for j in site_list:
                    if j[0][0][0]==i[0] and j[1]==site_env:
                        j[0].append(i)
                        site_saved=True
                        break
                if not site_saved:
                    site_list.append([[i],site_env])
            else:
                site_list.append([[i],site_env])
        for i in range(len(site_list)):
            site_list[i]=site_list[i][0]
        return site_list
    
    def sub(self,origin_site,new_site):
        poscar_list=[]
        for i in self.sites(site=origin_site,nei=2):
            if i[0][0]==origin_site:
                p=self.copy(True)
                origin_index=p.element.index(origin_site)
                if p.number[origin_index]==1:
                    del p.element[origin_index]
                    del p.number[origin_index]
                else:
                    p.number[origin_index]-=1
                origin_position=p.position[p.label.index(i[0])]
                del p.position[p.label.index(i[0])]
                if not new_site=="vac":
                    if not new_site in p.element:
                        p.element.insert(0,new_site)
                        p.number.insert(0,1)
                        p.position.insert(0,origin_position)
                        p.move(0.5-origin_position)
                    else:
                        new_index=p.element.index(new_site)
                        p.number[new_index]+=1
                        p.position.insert(sum(p.number[:new_index]),origin_position)
                p.label=[(i,k+1) for i,j in zip(p.element,p.number) for k in range(j)]
                poscar_list.append((i[0],p))
        return poscar_list
    
    def slab(self,vac_length,multi=5):
        import numpy as np
        # z方向扩展2倍
        self.number=[multi*i for i in self.number]
        c_length=np.linalg.norm(self.axis[-1])
        ratio=(c_length*multi+vac_length)/c_length
        for i in range(len(self.axis[-1])):
            self.axis[-1][i]*=ratio
        new_pos=[]
        for i in range(len(self.position)):
            for k in range(multi):
                new_pos.append(np.array([(self.position[i][j]+k)/ratio if j==2 else self.position[i][j] for j in range(3)]))
        self.position=new_pos
        return None
    
    def slab2(self,vac_length):
        import numpy as np
        self.number=[2*i for i in self.number]
        c_length=np.linalg.norm(self.axis[-1])
        ratio=(c_length*2+vac_length)/c_length
        for i in range(len(self.axis[-1])):
            self.axis[-1][i]*=ratio
        move_list=[]
        for i in self.layer[-1]:
            if abs(i[2][0]-i[2][1])<0.01:
                move_list.append((i[0],i[1]))
                
        new_pos=[]
        for i,j in zip(self.position,self.label):
            new_pos.append(np.array([i[j]/ratio if j==2 else i[j] for j in range(3)]))
            if j in move_list:
                new_pos.append(np.array([1+(i[j]-1)/ratio if j==2 else i[j] for j in range(3)]))
            else:
                new_pos.append(np.array([(i[j]+1)/ratio if j==2 else i[j] for j in range(3)]))
        self.position=new_pos
        return None
    
    def expand(self,tm):
        # transform_matrix
        import numpy as np
        self.axis=np.dot(tm,self.axis)
        #self.label=[(i,k+1) for i,j in zip(self.element,self.number) for k in range(j)]
        for i in range(len(self.number)):
            self.number[i]=0
        new_positions=[]
        for i,j in zip(self.label,self.position):
            for x in range(-3,4):
                for y in range(-3,4):
                    for z in range(-3,4):
                        a_new_position=(j+np.array([x,y,z])).dot(np.linalg.inv(tm))
                        if max(a_new_position)<1 and min(a_new_position)>=0:
                            new_positions.append(a_new_position)
                            self.number[self.element.index(i[0])]+=1
        self.position=new_positions
        return None

class job:
    def __init__(self):
        pass
    
    def read_job(self,string):
        import numpy as np
        self.state="wait"
        self.bundle=[]
        detail_list=csv_split(string)
        template_dict={
            "vasp":["type","structure","wavecar","occ","spd1","spd2","chg",
                    "kmesh","relax","soc","hdft","nele","spin","ldau","queue","core","encut","gga","mix","path"],
            "cc":["type","structure","vib1","vib2","ini","fin","step","path"],
            "cc2":["type","structure","vib1","vib2","ini","fin","step","path"],
            "para":["type","structure","parameter","ini","fin","step","path"],
            "distort":["type","structure","type","path"],
            "distort_sphere":["type","structure","distort_list","path"],
            "sub":["type","structure","origin_site","new_site","path"],
            "matproj":["type","formula","prop","path"],
            "atom":["type","element","path"],
            "slab":["type","structure","vac_length","path"],
            "slab2":["type","structure","vac_length","path"],
            "expand":["type","structure","tm","path"],
            "plot_split":["type","exm1","exm1_soc","grd_soc","spd","path"],
            "plot_split_cc":["type","cc_path","spd","path"],
            "plot_split_cc2":["type","cc1_path","cc2_path","spd","path"],
            "plot_ctl":["type","unitcell","pattern","path"],
            "plot_ctl3":["type","p1","m1","unitcell","path"],
            "ctl4":["type","system","level"],
            "plot_nei":["type","atom","range","path"],
            "plot_cc":["type","grd","ex","ex_at_grd","grd_at_ex","path"],
            "plot_cc2":["type","type","path"],
            "plot_pot":["type","pot","path"],
            "plot_band":["type","structure","path"],
            "show_procar":["type","structure","origin_structure","spd1","spd2"],
            "primcell":["type","structure","path"],
            "formation":["type","defect_list","host"]}
        self.detail={}
        for i,j in zip(template_dict[detail_list[0]],detail_list):
            if j=="FALSE":
                self.detail[i]=False
            elif j=="TRUE":
                self.detail[i]=True
            elif i in ["ini","fin","spin","encut","vac_length"]:
                self.detail[i]=float(j)
            elif i in ["core","range","step","relax"]:
                self.detail[i]=int(j)
            elif i=="nele":
                self.detail[i]=eval(j)
            elif i=="tm":
                self.detail[i]=[]
                j=j.split("_")
                for k in range(len(j)):
                    if not k%3:
                        self.detail[i].append([])
                    self.detail[i][-1].append(float(j[k]))
                self.detail[i]=np.array(self.detail[i])
            elif i=="ldau":
                self.detail[i]={}
                j=j.split("_")
                for k in range(len(j)):
                    if not k%4:
                        self.detail[i][j[k]]=[]
                        ok=k
                    else:
                        self.detail[i][j[ok]].append(j[k])
            else:
                self.detail[i]=j
        return None
    
    def copy(self):
        ajob=job()
        ajob.state=self.state
        ajob.bundle=self.bundle
        ajob.detail={}
        for i in self.detail:
            ajob.detail[i]=self.detail[i]
        return ajob

    
    def mk_kpoints_and_potcar(self):
        import re
        if self.detail["kmesh"]=="auto":
            linux_command("(echo 102;echo 2;echo 0.05)|vaspkit")
        elif self.detail["kmesh"]=="gamma":
            linux_command("(echo 102;echo 2;echo 0)|vaspkit")
        elif self.detail["kmesh"]=="line":
            linux_command("vaspkit -task 303")
            linux_command("cp KPATH.in KPOINTS")
        num_pattern=re.compile("([\d\.]+)")
        match_list=[]
        for i in linux_command("grep 'parameters from PSCTR are:' POTCAR -B1"):
            match_result=num_pattern.search(i)
            if match_result:
                match_list.append(float(match_result.group(1)))
        return match_list
    
    def get_encut(self):
        encut_list=[]
        with open("POTCAR","r") as f:
            for i in f.readlines():
                if "ENMAX" in i:
                    encut_list.append(float(i.split("=")[1].split(";")[0]))
        return max(encut_list)
            
    def mk_incar(self):
        import os
        import numpy as np
        cwd=os.getcwd()
        os.chdir(self.detail["path"])
        p=poscar()
        p.read("POSCAR")
        nel=self.detail["nele"]+np.dot(p.number,self.mk_kpoints_and_potcar())
        incar={
            "NCORE":int(np.power(self.detail["core"],0.5)/2)*2,
            "NELECT":nel,
            "LREAL":"Auto",
            "ENCUT":self.get_encut() if self.detail["encut"]<100 else self.detail["encut"],
            "ALGO":"All",
            "GGA":self.detail["gga"],
            "LORBIT":11,
            "NSW":0,
            "IBRION":-1,
            "EDIFF":1e-4,
            "NELM":50,
            "NELMIN":5,
            "ISMEAR":0,
            "SIGMA":0.05
            }
        if not self.detail["wavecar"]:
            incar["NELMDL"]=-5
        if self.detail["relax"]:
            if type(self.detail["relax"])==bool:
                incar["NSW"]=25
            else:
                incar["NSW"]=self.detail["relax"]
            # if "ez" in self.detail["path"] or "cz" in self.detail["path"] or "temp" in self.detail["path"]:
            #     incar["IBRION"]=1
            # else:
            #     incar["IBRION"]=2
            incar["IBRION"]=2
            # 计算介电张量
            if incar["NSW"]==1:
                incar["IBRION"]=8
                incar["LEPSILON"]=".TRUE."
                incar["NCORE"]=1
            incar["POTIM"]=0.5
            incar["ISIF"]=2
            incar["NELM"]=30
            incar["EDIFFG"]=-0.01
            incar["ALGO"]="Fast"
        if self.detail["soc"]:
            incar["LSORBIT"]="True"
            incar["SAXIS"]="0 0 1"
            # if "unitcell" in self.detail["path"]:
            #     incar["MAGMOM"]="600*0"
            # else:
            #     incar["MAGMOM"]="0 0 3 600*0"
        if self.detail["hdft"]:
            incar["ALGO"]="All"
            incar["LHFCALC"]=".TRUE."
            incar["PRECFOCK"]="Fast"
            if "pbe0" in self.detail["hdft"]:
                incar["AEXX"]=float(self.detail["hdft"].replace("pbe0",""))/100
            else:
                incar["AEXX"]=0.25
                if self.detail["hdft"]=="hse03":
                    incar["HFSCREEN"]=0.3
                elif self.detail["hdft"]=="hse06":
                    incar["HFSCREEN"]=0.2
        if self.detail["chg"]:
            incar["ICHARG"]=11
        os.chdir(cwd)
        
        if "slab" in self.detail["path"]:
            incar["LVHAR"]=".TRUE."
        
        if self.detail["occ"]:
            incar["ISMEAR"]=-2
            if self.detail["occ"]=="easy":
                incar["FERWE"]=str(int(nel)-1)+"*1 0 1 1000*0"
            else:
                fer=self.get_occ(0.5,0.5) if "half" in self.detail["path"] else self.get_occ()
                if not type(fer)==bool:
                    if len(fer)>1:
                        incar["FERWE"],incar["FERDO"]=fer
                    else:
                        incar["FERWE"]=fer[0]
                else:
                    if not fer:
                        incar["FERWE"]="failed"
                    
        os.chdir(self.detail["path"])
        if self.detail["mix"]:
            if self.detail["soc"]:
                incar["AMIX"]=0.2
                incar["BMIX"]=1e-5
                incar["AMIX_MAG"]=0.8
                incar["BMIX_MAG"]=1e-5
            else:
                incar["AMIX"]=0.2
                incar["BMIX"]=1e-5
        if not type(self.detail["spin"])==bool:
            incar["ISPIN"]=2
            incar["NUPDOWN"]=self.detail["spin"]
        elif self.detail["spin"]:
            incar["ISPIN"]=2
        if self.detail["ldau"]:
            incar["LDAU"]=".TRUE."
            incar["LDAUL"]=""
            incar["LDAUU"]=""
            incar["LDAUJ"]=""
            for i in p.element:
                if i in self.detail["ldau"]:
                    for j,k in zip(self.detail["ldau"][i],["LDAUL","LDAUU","LDAUJ"]):
                        incar[k]+=j+" "
                else:
                    for j,k in zip(["-1","0","0"],["LDAUL","LDAUU","LDAUJ"]):
                        incar[k]+=j+" "
        with open("INCAR","w") as f:
            for i in incar:
                f.write(i+"="+str(incar[i])+"\n")
        os.chdir(cwd)
        return None
    
    def runvasp(self,tag="normal"):
        import os
        import os.path
        cwd=os.getcwd()
        
        if tag=="rerun":
            pass
        elif tag=="normal":
            linux_command("mkdir "+self.detail["path"])
            linux_command("cp "+self.detail["structure"]+"/CONTCAR "+self.detail["path"]+"/POSCAR")
            if self.detail["wavecar"]:
                linux_command("cp "+self.detail["wavecar"]+"/WAVECAR "+self.detail["path"]+"/WAVECAR")
            if self.detail["chg"]:
                linux_command("cp "+self.detail["chg"]+"/CHGCAR "+self.detail["path"]+"/CHGCAR")
        # 统一节点设置
        # if self.detail["hdft"]:
        #     self.detail["queue"]="ckduan"
        #     self.detail["core"]=24
        # else:
        #     self.detail["queue"]="knl64"
        #     self.detail["core"]=64    
        self.mk_incar()
        os.chdir(self.detail["path"])
        hostname=linux_command("hostname")[0].replace("\n","")
        if hostname=="tc4600v3":
            if self.detail["soc"]:
                version="soc"
            elif self.detail["kmesh"]=="gamma":
                version="gamma"
            else:
                version="std"
            version_dir={"soc":"mpijob /opt/vasp/5.4.4/bin/vasp_ncl",
                          "gamma":"mpijob /opt/vasp/5.4.4/bin/vasp_gam",
                          "std":"mpijob /opt/vasp/5.4.4/bin/vasp_std"}
            
            if self.detail["queue"]=="ckduan":
                #白天模式
                #self.detail["queue"]+=" -m \"node197 node199 node200\""
                #夜间模式
                pass
            elif self.detail["queue"]=="test":
                self.detail["queue"]+=" -m \"k802 k804\""
            elif self.detail["queue"]=="knl64":
                self.detail["queue"]+=" -m \"knl3\""
            job_info=linux_command("bsub -q "+self.detail["queue"]+" -n "+str(self.detail["core"])+" -o %J.log "+version_dir[version])[0]
            
            print(job_info,end="")
            
            with open("jobid","a+") as q:
                q.write(job_info)
            with open("/home/phys/"+linux_command("whoami")[0].replace("\n","")+"/.jobs","a+") as q:
                q.write(job_info)
                # q.write(time.asctime(time.localtime(time.time()))+"\n")
                q.write(os.getcwd()+"\n")
                
        elif hostname=="login01":
            blist=["#!/bin/bash",
                   "#SBATCH -J test",
                   "#SBATCH -p compute",
                   "#SBATCH -n 48",
                   "#SBATCH -N 1",
                   "#SBATCH --time=24:00:00",
                   "#SBATCH -o %j.log",
                   "#SBATCH --qos=cpu",
                   "mpirun /share/softwares/vasp/from_SCC/"]
            if self.detail["soc"]:
                blist[-1]+="vasp_ncl"
            elif self.detail["kmesh"]=="gamma":
                blist[-1]+="vasp_gam"
            else:
                blist[-1]+="vasp_std"
                
            if os.path.isfile("vasp.sh"):
                linux_command("rm vasp.sh")
            
            with open("vasp.sh","w") as q:
                for i in blist:
                    q.write(i+"\n")
                    
            job_info=linux_command("sbatch vasp.sh")[0]
            
            print(job_info,end="")
            
            
            with open("jobid","a+") as q:
                q.write(job_info)
            with open("/share/home/ckduan/.jobs","a+") as q:
                q.write(job_info)
                q.write(os.getcwd()+"\n")
                    
        os.chdir(cwd)
        
        return None
    
    def get_occ(self,num1=0,num2=1):
        if type(self.detail["occ"])==bool and self.detail["occ"]:
            return True
        
        procar=read_procar(self.detail["occ"])
        
        if self.detail["spd1"]:
            spd1_list=[]
            # if "_" in self.detail["spd1"]:
            #     self.detail["spd1"]=tuple(self.detail["spd1"].split("_"))
            #     for j in procar[-1]:
            #         if self.detail["spd1"] in j["spd"]:
            #             spd1_list.append(j)
            # else:
            #     for j in procar[-1]:
            #         if j["spd"] and self.detail["spd1"] in list(zip(*j["spd"].keys()))[0]:
            #             spd1_list.append(j)
            
            if not "no_" in self.detail["spd1"]:
                spd1=tuple(self.detail["spd1"].split("_"))
                for j in procar[-1]:
                    if spd1 in j["spd"] and j["spd"][spd1]>0.05:
                        spd1_list.append(j)
            else:
                spd1=self.detail["spd1"].replace("no_","")
                for j in procar[-1]:
                    if not j["spd"]:
                        spd1_list.append(j)
                    else:
                        flag=True
                        for k in j["spd"]:
                            if spd1 in k and j["spd"][k]>0.05:
                                flag=False
                                break
                        if flag:
                            spd1_list.append(j)
                                
            for j in spd1_list[::-1]:
                if j["occ"]>=0.5:
                    procar[-1][procar[-1].index(j)]["occ"]=num1
                    break
            
        if self.detail["spd2"]:
            spd2_list=[]
            # if "_" in self.detail["spd2"]:
            #     self.detail["spd2"]=tuple(self.detail["spd2"].split("_"))
            #     for j in procar[0][::-1]:
            #         if self.detail["spd2"] in aaaaa["spd"]:
            #             spd2_list.append(j)
            # else:
            #     for j in procar[0][::-1]:
            #         if j["spd"] and self.detail["spd2"] in list(zip(*j["spd"].keys()))[0]:
            #             spd2_list.append(j)
            
            if not "no_" in self.detail["spd2"]:
                spd2=tuple(self.detail["spd2"].split("_"))
                for j in procar[0][::-1]:
                    if spd2 in j["spd"] and j["spd"][spd2]>0.1:
                        spd2_list.append(j)
            else:
                spd2=self.detail["spd2"].replace("no_","")
                for j in procar[0][::-1]:
                    if not j["spd"]:
                        spd2_list.append(j)
                    else:
                        flag=True
                        for k in j["spd"]:
                            if spd2 in k and j["spd"][k]>0.05:
                                flag=False
                                break
                        if flag:
                            spd2_list.append(j)
            
            for j in spd2_list[::-1]:
                if j["occ"]<=0.5:
                    procar[0][procar[0].index(j)]["occ"]=num2
                    break
        
        if len(procar)>1:
            occ=1
            count=0
            ferwe=""
            for j in range(len(procar[0])):
                if procar[0][j]["occ"]==occ:
                    count+=1
                else:
                    if count==1:
                        ferwe+=str(occ)+" "
                    else:
                        ferwe+=str(count)+"*"+str(occ)+" "
                    occ=procar[0][j]["occ"]
                    count=1
            ferwe+=str(count+1000)+"*"+str(occ)+" "
            
            occ=1
            count=0
            ferdo=""
            for j in range(len(procar[1])):
                if procar[1][j]["occ"]==occ:
                    count+=1
                else:
                    if count==1:
                        ferdo+=str(occ)+" "
                    else:
                        ferdo+=str(count)+"*"+str(occ)+" "
                    occ=procar[1][j]["occ"]
                    count=1
            ferdo+=str(count+1000)+"*"+str(occ)+" "
            
            occ_res=[ferwe,ferdo]
        else:
            occ=1
            count=0
            ferwe=""
            for j in range(len(procar[0])):
                if procar[0][j]["occ"]==occ:
                    count+=1
                else:
                    if count==1:
                        ferwe+=str(occ)+" "
                    else:
                        ferwe+=str(count)+"*"+str(occ)+" "
                    occ=procar[0][j]["occ"]
                    count=1
            ferwe+=str(count+1000)+"*"+str(occ)+" "
            occ_res=[ferwe]
        if self.detail["hdft"]:
            if "half" in self.detail["path"]:
                if occ_res[0].count("*")>3:
                    occ_res=False
            else:
                for j in occ_res:
                    if j.count("*")>2:
                        occ_res=False
                        break
        return occ_res

class jobs:
    def __init__(self):
        import os
        self.job_list=[]
        self.job_path_list=[]
        for i in os.listdir():
            if ".csv" in i:
                with open(i,"r",encoding='UTF-8-sig') as f:
                    for j in f.readlines():
                        #if not "#" in j and not "attempt" in j:
                        if not "#" in j:
                            ajob=job()
                            ajob.read_job(j)
                            # 去掉重复作业,不同策略除外
                            if not ajob.detail["type"]=="vasp":
                                if "path" in ajob.detail and ajob.detail["path"] in self.job_path_list:
                                    ajob=False
                            else:
                                for k in self.job_list:
                                    if k.detail==ajob.detail:
                                        ajob=False
                                        break
                            if ajob:
                                if "path" in ajob.detail:
                                    self.job_path_list.append(ajob.detail["path"])
                                self.job_list.append(ajob)
                                
    
    def find_job_by_path(self,path):
        i=False       
        for i in self.job_list:
            if "path" in i.detail:
                if i.detail["path"]==path:
                    return i
            
    def update_vasp(self,path):
        import re
        #import time
        import os
        
        ajob=self.find_job_by_path(path)
        listdir=os.listdir(path)
        # 获取最新的jobid
        if not "jobid" in os.listdir(path):
            log_list=[int(i.replace(".log","")) for i in os.listdir(path) if i.replace(".log","").isdigit()]
            if log_list:
                log_list.sort()
                for i in log_list:
                    linux_command("echo "+str(i)+" >> "+path+"/jobid")
        with open(path+"/jobid","r") as f:
            jobids=f.readlines()
        jobid=re.search("(\d+)",jobids[-1]).group(1)
        hostname=linux_command("hostname")[0].replace("\n","")
        if hostname=="tc4600v3":
            # 没有log文件
            if not jobid+".log" in listdir:
                # 没有OUTCAR
                if not "OUTCAR" in listdir:
                    state="run"
                # 作业超时
                # elif time.time()-os.path.getmtime(path+"/OUTCAR")>(3600 if "pbe0" in path else 600):
                #     linux_command("bkill "+jobid)
                #     state="fail"
                else:
                    state="run"
                    with open(path+"/OSZICAR","r") as f:
                        osz_lines=f.readlines()
                        #难以收敛
                        count=0
                        for i in osz_lines:
                            if ":  30    " in i:
                                count+=1
                            if count>5:
                                linux_command("bkill "+jobid)
                                return "fail"
                    
            else:
                with open(path+"/"+jobid+".log","r") as f:
                    log_lines=f.readlines()
                if "for stderr output of this job" in log_lines[-2]:
                    log_lines=log_lines[:-6]
                
                # Error EDDDAV: Call to ZHEGV failed. Returncode =   8 1   9
                if "Error EDDDAV: Call to ZHEGV failed." in log_lines[-1]:
                    if len(jobids)<5:
                        index=self.job_list.index(ajob)
                        self.job_list[index].detail["wavecar"]=False
                        self.job_list[index].detail["mix"]=True
                        linux_command("mv "+path+"/"+jobid+".log "+path+"/failed_"+jobid+".log")
                        linux_command("rm "+path+"/OUTCAR")
                        ajob.runvasp("rerun")
                        state="run"
                    else:
                        print("failcode1")
                        state="fail"
                
                if "VERY BAD NEWS!" in log_lines[-2]:
                    print("failcode2")
                    state="fail"
                        
                    
                # ZBRENT: fatal error in bracketing please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue
                # ZBRENT: fatal error: bracketing interval incorrect please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue
    
                elif "ZBRENT: fatal error" in log_lines[-3]:
                    if len(jobids)<5:
                        linux_command("cp "+path+"/CONTCAR "+path+"/POSCAR")
                        linux_command("mv "+path+"/"+jobid+".log "+path+"/failed_"+jobid+".log")
                        linux_command("rm "+path+"/OUTCAR")
                        ajob.runvasp("rerun")
                        state="run"     
                    else:
                        print("failcode3")
                        state="fail"
                
                elif "failed" in log_lines[-1] or "ERROR" in log_lines[-2] or "Error reading item 'FERWE' from file INCAR." in log_lines[-1]:
                    print("failcode4")
                    state="fail"
                
                else:
                    print(ajob.detail["relax"],not type(ajob.detail["relax"])==int,not "reached required accuracy - stopping structural energy minimisation" in log_lines[-2])
                    for i in log_lines:
                        if "# LSBATCH: User input" in i:
                            break
                        
                    # if "Exited with" in log_lines[log_lines.index(i)+4] or "TERM_" in log_lines[log_lines.index(i)+4]:
                    #     print("failcode5")
                    #     state="fail"
                    
                    if (ajob.detail["relax"] and not type(ajob.detail["relax"])==int and not "reached required accuracy - stopping structural energy minimisation" in log_lines[-2]):
                        if len(jobids)<5:
                            linux_command("cp "+path+"/CONTCAR "+path+"/POSCAR")
                            #linux_command("mv "+path+"/"+jobid+".log "+path+"/failed_"+jobid+".log")
                            linux_command("rm "+path+"/OUTCAR")
                            ajob.runvasp("rerun")
                            state="run"     
                        else:
                            print("failcode6")
                            state="fail"
                    else:
                        state="done"
        
        elif hostname=="login01":
            # 没有log文件
            if jobid in [[j for j in i.split(" ") if j][0] for i in linux_command("squeue")[1:]]:
                if "OSZICAR" in os.listdir(path):
                    with open(path+"/OSZICAR","r") as f:
                        osz_lines=f.readlines()
                        #难以收敛
                        count=0
                        for i in osz_lines:
                            if ":  30    " in i:
                                count+=1
                            if count>5:
                                linux_command("bkill "+jobid)
                                return "fail"
                state="run"
                
            # 作业超时
            # elif time.time()-os.path.getmtime(path+"/OUTCAR")>(3600 if "pbe0" in path else 600):
            #     linux_command("bkill "+jobid)
            #     state="fail"
            
            
            
            else:
                with open(path+"/"+jobid+".log","r") as f:
                    log_lines=f.readlines()
                
                # Error EDDDAV: Call to ZHEGV failed. Returncode =   8 1   9
                if "Error EDDDAV: Call to ZHEGV failed." in log_lines[-1] or "Error EDDDAV: Call to ZHEGV failed." in log_lines[-6]:
                    if len(jobids)<5:
                        index=self.job_list.index(ajob)
                        self.job_list[index].detail["wavecar"]=False
                        self.job_list[index].detail["mix"]=True
                        linux_command("mv "+path+"/"+jobid+".log "+path+"/failed_"+jobid+".log")
                        linux_command("rm "+path+"/OUTCAR")
                        ajob.runvasp("rerun")
                        state="run"
                    else:
                        state="fail"
                
                elif "VERY BAD NEWS!" in log_lines[-2]:
                    state="fail"
                
                elif "very serious problems" in log_lines[-2]:
                    state="fail"
                
                # ZBRENT: fatal error in bracketing please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue
                # ZBRENT: fatal error: bracketing interval incorrect please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue
    
                elif "ZBRENT: fatal error" in log_lines[-3] or "ZBRENT: fatal error" in log_lines[-8]:
                    if len(jobids)<5:
                        linux_command("cp "+path+"/CONTCAR "+path+"/POSCAR")
                        linux_command("mv "+path+"/"+jobid+".log "+path+"/failed_"+jobid+".log")
                        linux_command("rm "+path+"/OUTCAR")
                        ajob.runvasp("rerun")
                        state="run"     
                    else:
                        state="fail"
                
                elif "failed" in log_lines[-1] or "ERROR" in log_lines[-2] or "Error reading item 'FERWE' from file INCAR." in log_lines[-1]:
                    state="fail"
                
                elif ajob.detail["relax"] and not type(ajob.detail["relax"])==int:
                    
                    if "reached required accuracy - stopping structural energy minimisation" in log_lines[-2]:
                        state="done"
                    elif len(jobids)<5:
                        linux_command("cp "+path+"/CONTCAR "+path+"/POSCAR")
                        linux_command("mv "+path+"/"+jobid+".log "+path+"/failed_"+jobid+".log")
                        linux_command("rm "+path+"/OUTCAR")
                        ajob.runvasp("rerun")
                        state="run"
                    else:
                        state="fail"
                else:
                    state="done"
        return state
    
    def update(self,file,delay):
        import os
        import re
        import numpy as np
        jobid_template=re.compile("(\d+)")
        for i in self.job_list:
            print(i.detail)
            #print(i.bundle)
            #print(i.state)
            if i.detail["type"]=="vasp":
                
                # 判断vasp计算中的POSCAR对称性
                # if os.path.isdir(i.detail["path"]):
                #     pwd=os.getcwd()
                #     os.chdir(i.detail["path"])
                #     for j in linux_command("vaspkit -task 601"):
                #         if "Point Group" in j:
                #             file.write("point "+j.split("[")[1].split("]")[0]+" "+i.detail["path"]+"\n")
                #             break
                #     os.chdir(pwd)
                
                if i.state in ["wait","run"]:
                    # 更新vasp的bundle
                    job_structure=self.find_job_by_path(i.detail["structure"])
                    if i.detail["path"] in self.job_path_list:
                        # 不是子任务
                        if job_structure.bundle:
                            # bundle
                            if i.bundle:
                                state_list=[self.find_job_by_path(j).state for j in i.bundle]
                                if "fail" in state_list:
                                    i.state="fail"
                                elif set(state_list)=={"done"}:
                                    i.state="done"
                            elif os.path.isdir(i.detail["path"]):
                                # bundle丢失
                                for j in os.listdir(i.detail["structure"]):
                                    self.job_list.append(i.copy())
                                    self.job_list[-1].detail["structure"]=self.job_list[-1].detail["structure"]+"/"+j
                                    self.job_list[-1].detail["path"]=self.job_list[-1].detail["path"]+"/"+j
                                    #self.job_list[-1].state=self.update_vasp(self.job_list[-1].detail["path"])
                                    i.bundle.append(self.job_list[-1].detail["path"])
                                    i.state="run"
                                    
                            elif (job_structure.state=="done" 
                                  and (not i.detail["wavecar"] or self.find_job_by_path(i.detail["wavecar"]).state=="done") 
                                  and (type(i.detail["occ"])==bool or self.find_job_by_path(i.detail["occ"]).state=="done")
                                  and (not i.detail["chg"] or self.find_job_by_path(i.detail["chg"]).state=="done")):
                                # bundle未提交
                                linux_command("mkdir "+i.detail["path"])
                                for j in os.listdir(i.detail["structure"]):
                                    self.job_list.append(i.copy())
                                    self.job_list[-1].detail["structure"]=self.job_list[-1].detail["structure"]+"/"+j
                                    self.job_list[-1].detail["path"]=self.job_list[-1].detail["path"]+"/"+j
                                    self.job_list[-1].runvasp()
                                    self.job_list[-1].state="run"
                                    i.bundle.append(self.job_list[-1].detail["path"])
                                    i.state="run"
                        else:
                            # no bundle
                            if os.path.isdir(i.detail["path"]):
                                i.state=self.update_vasp(i.detail["path"])
                            elif (job_structure.state=="done" 
                                  and (not i.detail["wavecar"] or self.find_job_by_path(i.detail["wavecar"]).state=="done") 
                                  and (type(i.detail["occ"])==bool or i.detail["occ"]=="easy" or self.find_job_by_path(i.detail["occ"]).state=="done")
                                  and (not i.detail["chg"] or self.find_job_by_path(i.detail["chg"]).state=="done")):
                                i.runvasp()
                                i.state="run"
                    # 路径存在
                    elif os.path.isdir(i.detail["structure"]):
                        # 子任务
                        if os.path.isdir(i.detail["path"]):
                            i.state=self.update_vasp(i.detail["path"])
                        else:
                            i.runvasp()
                            i.state="run"
                    
            # 非vasp
            elif i.state=="wait":
                # 1.prepared?
                # 2.作业路径
                # 3.bundle?
                # 4.p.read
                # 5.done
                if i.detail["type"]=="matproj":
                    if os.path.isdir(i.detail["path"]):
                        if not "CONTCAR" in os.listdir(i.detail["path"]):
                            for j in os.listdir(i.detail["path"]):
                                i.bundle.append(i.detail["path"]+"/"+j)
                        i.state="done"
                    else:
                        linux_command("mkdir "+i.detail["path"])
                        i.bundle=matproj(i.detail["formula"],i.detail["prop"],i.detail["path"])
                        i.state="done"
                    
                elif i.detail["type"]=="atom":
                    if os.path.isdir(i.detail["path"]):
                        i.state="done"
                    else:
                        p=poscar()
                        p.atom(i.detail["element"])
                        linux_command("mkdir "+i.detail["path"])
                        p.write(i.detail["path"]+"/CONTCAR")
                        i.state="done"
                

                elif i.detail["type"]=="formation" and not delay:
                    chempot_sample(i.detail["formula"],sample_times=1000)
                    # Cl_rich={"Cl":-2.31253666,"Bi":-8.605815231,"Cs":-5.113613271,"Na":-5.136552897,"vac":0}
                    # Cl_moderate={"Cl":-3.05553666,"Bi":-6.204012487,"Cs":-4.485310508,"Na":-4.335106903,"vac":0}
                    # Cl_poor={"Cl":-3.91253666,"Bi":-3.894156521,"Cs":-3.48223995,"Na":-3.5064167,"vac":0}
                    
                    # Cl_rich={"Cl":-1.68453666,"Bi":-8.651844939,"Cs":-5.126672667,"Na":-5.143896817,"vac":0}
                    # Cl_moderate={"Cl":-2.35253666,"Bi":-6.662208715,"Cs":-4.476663706,"Na":-4.427074424,"vac":0}
                    # Cl_poor={"Cl":-3.04553666,"Bi":-4.714041505,"Cs":-3.71842408,"Na":-3.735319215,"vac":0}
                    
                    # cp=[['Cl_poor',
                    #      {'Cl': -3.19453666,
                    #       'Bi': -4.70088538,
                    #       'Cs': -3.734779862,
                    #       'Na': -3.594577483,
                    #       'vac': 0}],
                    #     ['Cl_moderate',
                    #      {'Cl': -2.41553666,
                    #       'Bi': -6.993367852,
                    #       'Cs': -4.528673198,
                    #       'Na': -4.392762595,
                    #       'vac': 0}],
                    #     ['Cl_rich',
                    #      {'Cl': -1.70553666,
                    #       'Bi': -9.040191304,
                    #       'Cs': -5.278674859,
                    #       'Na': -5.105603604,
                    #       'vac': 0}]]
                    
                    # cp=[['Cl_poor_ClNa_BiCl6Cs3',
                    #      {'Cl': -3.19453666,
                    #       'Bi': -4.78996531,
                    #       'Cs': -3.695402372,
                    #       'Na': -3.584525293,
                    #       'O':-4.947961005-2.7987541192998857,
                    #       'vac': 0}],
                    #     ['Cl_poor_ClNa_Bi_rich',
                    #      {'Cl': -3.19453666,
                    #       'Bi': -4.69796531,
                    #       'Cs': -3.741402372,
                    #       'Na': -3.584525293,
                    #       'O':-4.947961005-2.7987541192998857,
                    #       'vac': 0}],
                    #     ['Cl_poor_BiCl6Cs3_Bi_rich',
                    #      {'Cl': -3.19453666,
                    #       'Bi': -4.69796531,
                    #       'Cs': -3.726402372,
                    #       'Na': -3.615525293,
                    #       'O':-4.947961005-2.7987541192998857,
                    #       'vac': 0}],
                    #     ['Cl_moderate_ClNa_BiCl6Cs3',
                    #      {'Cl': -2.41553666,
                    #       'Bi': -7.12996531,
                    #       'Cs': -4.476402372,
                    #       'Na': -4.372525293,
                    #       'O':-4.947961005-1.503403011461916,
                    #       'vac': 0}],
                    #     ['Cl_moderate_ClNa_Bi2Cl9Cs3',
                    #      {'Cl': -2.41553666,
                    #       'Bi': -6.77196531,
                    #       'Cs': -4.649402372,
                    #       'Na': -4.372525293,
                    #       'O':-4.947961005-1.503403011461916,
                    #       'vac': 0}],
                    #     ['Cl_moderate_BiCl6Cs3_Bi2Cl9Cs3',
                    #      {'Cl': -2.41553666,
                    #       'Bi': -6.94496531,
                    #       'Cs': -4.532402372,
                    #       'Na': -4.421525293,
                    #       'O':-4.947961005-1.503403011461916,
                    #       'vac': 0}],
                    #     ['Cl_rich_ClNa_BiCl6Cs3',
                    #      {'Cl': -1.70553666,
                    #       'Bi': -9.24796531,
                    #       'Cs': -5.181402372,
                    #       'Na': -5.075525293,
                    #       'O':-4.947961005-0.46878133412934936,
                    #       'vac': 0}],
                    #     ['Cl_rich_ClNa_Bi2Cl9Cs3',
                    #      {'Cl': -1.70553666,
                    #       'Bi': -8.89896531,
                    #       'Cs': -5.364402372,
                    #       'Na': -5.075525293,
                    #       'O':-4.947961005-0.46878133412934936,
                    #       'vac': 0}],
                    #     ['Cl_rich_BiCl6Cs3_Bi2Cl9Cs3',
                    #      {'Cl': -1.70553666,
                    #       'Bi': -9.08196531,
                    #       'Cs': -5.239402372,
                    #       'Na': -5.141525293,
                    #       'O':-4.947961005-0.46878133412934936,
                    #       'vac': 0}]]
                    
                                        # cp=[['Cl_poor_ClNa_BiCl6Cs3',
                    #      {'Cl': -3.19453666,
                    #       'Bi': -4.78996531,
                    #       'Cs': -3.695402372,
                    #       'Na': -3.584525293,
                    #       'O':-4.947961005-2.7987541192998857,
                    #       'vac': 0}],
                    #     ['Cl_poor_ClNa_Bi_rich',
                    #      {'Cl': -3.19453666,
                    #       'Bi': -4.69796531,
                    #       'Cs': -3.741402372,
                    #       'Na': -3.584525293,
                    #       'O':-4.947961005-2.7987541192998857,
                    #       'vac': 0}],
                    #     ['Cl_poor_BiCl6Cs3_Bi_rich',
                    #      {'Cl': -3.19453666,
                    #       'Bi': -4.69796531,
                    #       'Cs': -3.726402372,
                    #       'Na': -3.615525293,
                    #       'O':-4.947961005-2.7987541192998857,
                    #       'vac': 0}],
                    #     ['Cl_moderate_ClNa_BiCl6Cs3',
                    #      {'Cl': -2.41553666,
                    #       'Bi': -7.12996531,
                    #       'Cs': -4.476402372,
                    #       'Na': -4.372525293,
                    #       'O':-4.947961005-1.503403011461916,
                    #       'vac': 0}],
                    #     ['Cl_moderate_ClNa_Bi2Cl9Cs3',
                    #      {'Cl': -2.41553666,
                    #       'Bi': -6.77196531,
                    #       'Cs': -4.649402372,
                    #       'Na': -4.372525293,
                    #       'O':-4.947961005-1.503403011461916,
                    #       'vac': 0}],
                    #     ['Cl_moderate_BiCl6Cs3_Bi2Cl9Cs3',
                    #      {'Cl': -2.41553666,
                    #       'Bi': -6.94496531,
                    #       'Cs': -4.532402372,
                    #       'Na': -4.421525293,
                    #       'O':-4.947961005-1.503403011461916,
                    #       'vac': 0}],
                    #     ['Cl_rich_ClNa_BiCl6Cs3',
                    #      {'Cl': -1.70553666,
                    #       'Bi': -9.24796531,
                    #       'Cs': -5.181402372,
                    #       'Na': -5.075525293,
                    #       'O':-4.947961005-0.46878133412934936,
                    #       'vac': 0}],
                    #     ['Cl_rich_ClNa_Bi2Cl9Cs3',
                    #      {'Cl': -1.70553666,
                    #       'Bi': -8.89896531,
                    #       'Cs': -5.364402372,
                    #       'Na': -5.075525293,
                    #       'O':-4.947961005-0.46878133412934936,
                    #       'vac': 0}],
                    #     ['Cl_rich_BiCl6Cs3_Bi2Cl9Cs3',
                    #      {'Cl': -1.70553666,
                    #       'Bi': -9.08196531,
                    #       'Cs': -5.239402372,
                    #       'Na': -5.141525293,
                    #       'O':-4.947961005-0.46878133412934936,
                    #       'vac': 0}]]
                    # cp=[['Sn',{'Cl': -1.84853666-0.00022831100254822696,
                    #            'Cs': -0.8954023720689656-5.031111461912579,
                    #            'Sn': -4.009571855-7.237553751021622,
                    #            'vac': 0}],
                    #     ['Sn',{'Cl': -1.84853666-0.00022831100254822696,
                    #            'Cs': -0.8954023720689656-5.031111461912579,
                    #            'Sn': -4.009571855-7.237553751021622,
                    #            'vac': 0}],
                    #     ['Sn',{'Cl': -1.84853666-0.00022831100254822696,
                    #            'Cs': -0.8954023720689656-5.031111461912579,
                    #            'Sn': -4.009571855-7.237553751021622,
                    #            'vac': 0}],
                    #     ['Sn',{'Cl': -1.84853666-0.00022831100254822696,
                    #            'Cs': -0.8954023720689656-5.031111461912579,
                    #            'Sn': -4.009571855-7.237553751021622,
                    #            'vac': 0}],
                    #     ['Sn',{'Cl': -1.84853666-0.00022831100254822696,
                    #            'Cs': -0.8954023720689656-5.031111461912579,
                    #            'Sn': -4.009571855-7.237553751021622,
                    #            'vac': 0}],
                    #     ['Sn',{'Cl': -1.84853666-0.00022831100254822696,
                    #            'Cs': -0.8954023720689656-5.031111461912579,
                    #            'Sn': -4.009571855-7.237553751021622,
                    #            'vac': 0}],]
                    cp1={(-6.557302319011928e-06, -4.829699517449182, -7.6417081621497935): ['Cl_max'], (-1.862286378977078, -2.6565572883836097, -0.8143136902323818): ['Cl_min'], (-1.861918334199625, -2.6483476389742884, -0.8329412577157427): ['Cs_max'], (-0.0004413290205943794, -5.08174127271423, -7.1350160213100455): ['Cs_min'], (-1.7575838548057412, -2.97086307264433, -0.8139172667389629): ['Sn_max'], (-0.0020041762071920165, -4.5085479430552216, -8.272025597508476): ['Sn_min']}
                    cp=[]
                    for j in cp1:
                        cp.append([cp1[j][0],{"Cl":-1.23453666+j[0],
                                              "Cs":-0.8954023720689656+j[1],
                                              "Sn":-4.009571855+j[2],
                                              'vac': 0}])
                        
                        
                    # for j in cp:
                    #     j[1]["O"]=-4.947961005
                    
                    pot={"Cs":0.388,"Na":0.142,"Bi":0.303,"Cl":0.15,"vac":0,"O":0.059,"Sn":0.294}
                    # for j in pot:
                    #     pot[j]/=4
                    
                    for k in cp:
                        defect_ctl={}
                        #j=self.find_job_by_path("0Cs2NaBiCl6/relax_grd_soc")
                        j=self.find_job_by_path("0Cs2SnCl6/relax_unitcell")
                        grd_energy=energy(j.detail["path"])-charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                        for j in self.job_list:
                            #if j.detail["type"]=="vasp" and not "grd" in j.detail["path"]:
                            if j.detail["type"]=="vasp" and not "unitcell" in j.detail["path"]:
                                if not defect_ctl or not j.detail["path"].split("_")[1] in defect_ctl:
                                    defect_ctl[j.detail["path"].split("_")[1]]=[]
                                asub=self.find_job_by_path(j.detail["structure"])
                                # "system":j.detail["path"].split("_")[1],
                                # "nele":-j.detail["nele"],
                                # "energy":energy(j.detail["path"])
                                #                         -charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                                #                         -grd_energy
                                #                         +Cl_rich[asub.detail["origin_site"]]
                                #                         -Cl_rich[asub.detail["new_site"]]
                                #                         -j.detail["nele"]*(-pot[asub.detail["origin_site"]]+pot[asub.detail["new_site"]])),
                                #                         j.detail["path"].split("_")[1]])
                                defect_ctl[j.detail["path"].split("_")[1]].append([-j.detail["nele"],
                                                        (energy(j.detail["path"])
                                                        -charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                                                        -grd_energy
                                                        +(k[1][asub.detail["origin_site"]] 
                                                          if not "-" in asub.detail["origin_site"] 
                                                          else sum([k[1][m] for m in asub.detail["origin_site"].split("-")]))
                                                        -(k[1][asub.detail["new_site"]]
                                                          if not "-" in asub.detail["new_site"] 
                                                          else sum([k[1][m] for m in asub.detail["new_site"].split("-")]))
                                                        -j.detail["nele"]*(-(pot[asub.detail["origin_site"]] 
                                                          if not "-" in asub.detail["origin_site"] 
                                                          else sum([pot[m] for m in asub.detail["origin_site"].split("-")]))
                                                                           +(pot[asub.detail["new_site"]]
                                                          if not "-" in asub.detail["new_site"] 
                                                          else sum([pot[m] for m in asub.detail["new_site"].split("-")])))),
                                                        j.detail["path"].split("_")[1]])
                        
                        #plot_ctl5([0.26434953,3.29911075],list(defect_ctl.values()),0.05,k[0])
                        plot_ctl5([-0.58145701,2.06722122],list(defect_ctl.values()),0.05,k[0])
                    
                    # defect_ctl={}
                    # j=self.find_job_by_path("0Cs2NaBiCl6/relax_grd_soc")
                    # grd_energy=energy(j.detail["path"])-charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                    # for j in self.job_list:
                    #     if j.detail["type"]=="vasp" and not "grd" in j.detail["path"]:
                    #         if not defect_ctl or not j.detail["path"].split("_")[1] in defect_ctl:
                    #             defect_ctl[j.detail["path"].split("_")[1]]=[]
                    #         asub=self.find_job_by_path(j.detail["structure"])
                    #         defect_ctl[j.detail["path"].split("_")[1]].append([-j.detail["nele"],
                    #                                 (energy(j.detail["path"])
                    #                                 -charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                    #                                 -grd_energy
                    #                                 +Cl_poor[asub.detail["origin_site"]]
                    #                                 -Cl_poor[asub.detail["new_site"]]
                    #                                 -j.detail["nele"]*(-pot[asub.detail["origin_site"]]+pot[asub.detail["new_site"]])),
                    #                                 j.detail["path"].split("_")[1]])
                    
                    # plot_ctl4([0.26434953,3.29911075],list(defect_ctl.values()),"Cl_poor")
                    
                    # defect_ctl={}
                    # j=self.find_job_by_path("0Cs2NaBiCl6/relax_grd_soc")
                    # grd_energy=energy(j.detail["path"])-charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                    # for j in self.job_list:
                    #     if j.detail["type"]=="vasp" and not "grd" in j.detail["path"]:
                    #         if not defect_ctl or not j.detail["path"].split("_")[1] in defect_ctl:
                    #             defect_ctl[j.detail["path"].split("_")[1]]=[]
                    #         asub=self.find_job_by_path(j.detail["structure"])
                    #         defect_ctl[j.detail["path"].split("_")[1]].append([-j.detail["nele"],
                    #                                 (energy(j.detail["path"])
                    #                                 -charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                    #                                 -grd_energy
                    #                                 +Cl_moderate[asub.detail["origin_site"]]
                    #                                 -Cl_moderate[asub.detail["new_site"]]
                    #                                 -j.detail["nele"]*(-pot[asub.detail["origin_site"]]+pot[asub.detail["new_site"]])),
                    #                                 j.detail["path"].split("_")[1]])
                    
                    # plot_ctl4([0.26434953,3.29911075],list(defect_ctl.values()),"Cl_moderate")
                    
                    # defect_ctl={}
                    # j=self.find_job_by_path("0Cs2NaBiCl6/relax_grd_soc")
                    # grd_energy=energy(j.detail["path"])-charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                    # for j in self.job_list:
                    #     if j.detail["type"]=="vasp" and not "grd" in j.detail["path"]:
                    #         if not defect_ctl or not j.detail["path"].split("_")[1] in defect_ctl:
                    #             defect_ctl[j.detail["path"].split("_")[1]]=[]
                    #         asub=self.find_job_by_path(j.detail["structure"])
                    #         defect_ctl[j.detail["path"].split("_")[1]].append([-j.detail["nele"],
                    #                                 (energy(j.detail["path"])
                    #                                 -charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                    #                                 -grd_energy
                    #                                 +Cl_moderate2[asub.detail["origin_site"]]
                    #                                 -Cl_moderate2[asub.detail["new_site"]]
                    #                                 -j.detail["nele"]*(-pot[asub.detail["origin_site"]]+pot[asub.detail["new_site"]])),
                    #                                 j.detail["path"].split("_")[1]])
                    
                    # plot_ctl4([0.26434953,3.29911075],list(defect_ctl.values()),"Cl_moderate2")
                
                # elif i.detail["type"]=="formation" and not delay:
                #     a=[[{"Cl":-1.938392256,"Cs":	-5.874532685,"Sn":	-10.81332419,"vac":0},"Sn rich    "],
                #      [{"Cl":-2.020088584,"Cs":	-5.241919459,"Sn":	-11.58837268,"vac":0},"Cs rich    "],
                #      [{"Cl":-3.636359261,"Cs":	-3.665097587,"Sn":	-5.044392358,"vac":0},"Cl Cs rich "],
                #      [{"Cl":-3.358611001,"Cs":	-4.414505915,"Sn":	-5.212065262,"vac":0},"Cl Sn rich "]]
                #     pot={"Cs":0.388,"Sn":0.294,"Cl":0.156,"vac":0}
                    
                #     for z in a:
                #         defect_ctl={}
                #         j=self.find_job_by_path("0Cs2SnCl6/relax_unitcell")
                #         grd_energy=energy(j.detail["path"])-charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                #         print
                #         for j in self.job_list:
                #             if j.detail["type"]=="vasp" and not "unitcell" in j.detail["path"]:
                #                 if not defect_ctl or not j.detail["path"].split("_")[1] in defect_ctl:
                #                     defect_ctl[j.detail["path"].split("_")[1]]=[]
                #                 asub=self.find_job_by_path(j.detail["structure"])
                #                 defect_ctl[j.detail["path"].split("_")[1]].append([-j.detail["nele"],
                #                                        (energy(j.detail["path"])
                #                                         -charge_corr_dict[j.detail["path"].split("/")[0]]*j.detail["nele"]**2
                #                                         -grd_energy
                #                                         +z[0][asub.detail["origin_site"]]
                #                                         -z[0][asub.detail["new_site"]]
                #                                         -j.detail["nele"]*(-pot[asub.detail["origin_site"]]+pot[asub.detail["new_site"]])),
                #                                        j.detail["path"].split("_")[1]])
                        
                #         plot_ctl4([-0.58145701,2.06722122],list(defect_ctl.values()),z[-1])
                    
                            
                            
                            
                            
                
                elif i.detail["type"]=="show_procar":
                    if (self.find_job_by_path(i.detail["structure"]).state=="done" 
                        and self.find_job_by_path(i.detail["origin_structure"]).state=="done"):
                        
                        job_structure=self.find_job_by_path(i.detail["structure"])
                        file.write(i.detail["structure"]+"\n")
                        # 和get_occ 相同
                        procar=read_procar(i.detail["structure"])
                        start=False
                        end=False
                        if i.detail["spd1"]:
                            if "_" in i.detail["spd1"]:
                                i.detail["spd1"]=tuple(i.detail["spd1"].split("_"))
                            for j in procar[-1]:
                                if j["occ"]<=0.5:
                                    # if type(i.detail["spd1"])==tuple:
                                    #     if not i.detail["spd1"] in j["spd"]:
                                    #         job_structure.state="wrong"
                                    # else:
                                    #     if not j["spd"] or not i.detail["spd1"] in list(zip(*j["spd"].keys()))[0]:
                                    #         job_structure.state="wrong"
                                    # start=procar[-1].index(j)
                                    
                                    if not "no" in i.detail["spd1"]:
                                        if not i.detail["spd1"] in j["spd"]:
                                            job_structure.state="wrong"
                                    else:
                                        if j["spd"]:
                                            for k in j["spd"]:
                                                if i.detail["spd1"][1] in k and j["spd"][k]>0.05:
                                                    job_structure.state="wrong"
                                                    break
                                    
                                    start=procar[-1].index(j)
                                    break
                        
                        if i.detail["spd2"]:
                            if "_" in i.detail["spd2"]:
                                i.detail["spd2"]=tuple(i.detail["spd2"].split("_"))
                            for j in procar[0][::-1]:
                                if j["occ"]>=0.5:
                                    # if type(i.detail["spd2"])==tuple:
                                    #     if not i.detail["spd2"] in j["spd"]:
                                    #         job_structure.state="wrong"
                                    # else:
                                    #     if not j["spd"] or not i.detail["spd2"] in list(zip(*j["spd"].keys()))[0]:
                                    #         job_structure.state="wrong"
                                    
                                    if not "no" in i.detail["spd2"]:
                                        if not i.detail["spd2"] in j["spd"]:
                                            job_structure.state="wrong"
                                    else:
                                        if j["spd"]:
                                            for k in j["spd"]:
                                                if i.detail["spd2"][1] in k and j["spd"][k]>0.05:
                                                    job_structure.state="wrong"
                                                    break
                                            
                                    end=procar[0].index(j)
                                    break
                        
                        real_start=start-5 if start else end-5
                        end=end+5 if end else start+5
                        
                        # 检查结构
                        distort_poscar=poscar()
                        origin_poscar=poscar()
                        origin_path=i.detail["origin_structure"]
                            
                        # 根据预测的畸变计算控制占据
                        if not "attempt" in i.detail["structure"] and "relax" in i.detail["structure"]:
                            if i.detail["structure"].replace("relax_","")+"_attempt" in self.job_path_list:
                                if (job_structure.state in ["wrong","fail"] or 
                                    self.find_job_by_path(job_structure.detail["structure"]).state in ["wrong","fail"]):
                                    ajob=self.find_job_by_path(i.detail["structure"].replace("relax_","")+"_attempt")
                                    ajob.detail["structure"]=ajob.detail["structure"].replace("_attempt","")
                                elif job_structure.state in ["done","run","wait"]:
                                    for j in self.job_path_list:
                                        if (i.detail["structure"].replace("relax_","")+"_attempt" in j or
                                            i.detail["structure"]+"_attempt" in j):
                                            if os.path.isdir(j):
                                                if "jobid" in os.listdir(j):
                                                    with open(j+"/jobid","r") as f:
                                                        linux_command("bkill "+jobid_template.search(f.readlines()[-1]).group(1))
                                                linux_command("mv "+j+" "+j+"_nomore")
                                                self.find_job_by_path(j).state="wait"
                                                
                        origin_poscar.read(origin_path+"/CONTCAR")
                        distort_poscar.read(i.detail["structure"]+"/CONTCAR")
                        distortion=distort_poscar-origin_poscar
                        # 对于低对称情况缺陷也可能移动,故不能直接使用poscar.distance()
                        nearest_distortion=np.linalg.norm(
                            distortion.position[
                                distort_poscar.label.index(
                                    distort_poscar.distance(
                                        distort_poscar.label[0],False)[1][0])].dot(distortion.axis))
                        
                        distort_list=[]
                        current_element=False
                        neis=0
                        ele_copy=False
                        if "unitcell" in i.detail["structure"]:
                            ban_list=["Br","F","Cl","Cs","Na","K","Ag"]
                            ele_copy=distort_poscar.element[:]
                            for j in ban_list:
                                if len(ele_copy)==1:
                                    break
                                elif j in ele_copy:
                                    ele_copy.remove(j)
                        ele_copy=(ele_copy[0],1) if ele_copy else distort_poscar.label[0]
                        for j in distort_poscar.distance(ele_copy,False):
                            if not current_element==j[0][0]:
                                neis+=1
                                if neis>3:
                                    break
                                else:
                                    current_element=j[0][0]
                            temp_distort_list=[]
                            file.write("%-12s" %str(j[0]))
                            temp_distort_list.append(j[0][0])
                            abs_origin=(origin_poscar.position[distort_poscar.label.index(j[0])]
                                    -origin_poscar.position[distort_poscar.label.index(ele_copy)])
                            for k in range(len(abs_origin)):
                                if abs_origin[k]>0.5:
                                    abs_origin[k]-=1
                                elif abs_origin[k]<=-0.5:
                                    abs_origin[k]+=1
                            abs_origin=abs_origin.dot(origin_poscar.axis)
                            for k in abs_origin:
                                file.write("%-8s" %str(round(k,3)))
                                
                            file.write("%-8s" %str(round(np.linalg.norm(abs_origin),3)))
                            temp_distort_list.append(round(np.linalg.norm(abs_origin),3))
                            abs_distortion=distortion.position[distort_poscar.label.index(j[0])].dot(distortion.axis)
                            for k in abs_distortion:
                                file.write("%-8s" %str(round(k/nearest_distortion,1)))
                                
                            file.write("%-8s" %str(round(np.linalg.norm(abs_distortion),3)))
                            
                            if np.linalg.norm(abs_origin):
                                file.write("%-8s" %str(round(abs_origin.dot(abs_distortion)/np.linalg.norm(abs_origin),3)))
                                temp_distort_list.append(round(abs_origin.dot(abs_distortion)/np.linalg.norm(abs_origin),3))
                            else:
                                file.write("%-8s" %str(round(np.linalg.norm(abs_distortion),3)))
                            #if not temp_distort_list in distort_list:
                            distort_list.append(temp_distort_list)
                            file.write("\n")
                        del distort_list[0]
                        file.write(job_structure.state+" distort "+i.detail["structure"]+" ")
                        if distort_list[0]==distort_list[5]:
                            distort_mode="A"
                        elif distort_list[0]==distort_list[3] and distort_list[4]==distort_list[5]:
                            distort_mode="Ee"
                        elif distort_list[0]==distort_list[1] and distort_list[2]==distort_list[5]:
                            distort_mode="Ec"
                        else:
                            distort_mode="unk"
                        file.write(distort_mode+" ")
                        prev=False
                        for j in distort_list:
                            if not prev==j:
                                for k in j:
                                    file.write(str(k)+" ")
                                prev=j
                        file.write("\n")
                    
                        # attempt成功则替换
                        if ("attempt" in i.detail["structure"]
                            and "relax" in i.detail["structure"]
                            and job_structure.state=="done"
                            and self.find_job_by_path(i.detail["structure"].split("_attempt")[0]).state in ["fail","wrong"]):
                            linux_command("mv "+i.detail["structure"].split("_attempt")[0]
                                  +" "+i.detail["structure"].split("_attempt")[0]+"_fail")
                            linux_command("mv "+i.detail["structure"]+" "+i.detail["structure"].split("_attempt")[0])
                            self.find_job_by_path(i.detail["structure"].split("_attempt")[0]).state="done"

                        # 检查激发态
                        for j in range(end)[real_start:]:
                            for k in procar:
                                file.write(str(k[j])+"\n")
                        
                        i.state="done"
                        
                    elif ("relax" in i.detail["structure"] 
                          and not "attempt" in i.detail["structure"] 
                          and self.find_job_by_path(i.detail["structure"]).state=="fail"
                          and i.detail["structure"].replace("relax_","")+"_attempt" in self.job_path_list):
                        ajob=self.find_job_by_path(i.detail["structure"].replace("relax_","")+"_attempt")
                        ajob.detail["structure"]=ajob.detail["structure"].replace("_attempt","")
                        
                        i.state="done"
                        
                
                elif i.detail["type"]=="sub":
                    if os.path.isdir(i.detail["path"]):
                        listdir=os.listdir(i.detail["path"])
                        if not "CONTCAR" in listdir:
                            for j in listdir:
                                i.bundle.append(i.detail["path"]+"/"+j)
                        i.state="done"
                        
                    elif self.find_job_by_path(i.detail["structure"]).state=="done":
                        p=poscar()
                        if self.find_job_by_path(i.detail["structure"]).bundle:
                            sub_list=[]
                            for j in self.find_job_by_path(i.detail["structure"]).bundle:
                                p.read(j+"/CONTCAR")
                                sub_list+=p.sub(i.detail["origin_site"],i.detail["new_site"])
                            if not os.path.isdir(i.detail["path"]):
                                linux_command("mkdir "+i.detail["path"])
                            if len(sub_list)==1:
                                sub_list[0][1].write(i.detail["path"]+"/CONTCAR")
                            else:
                                file.write("\n"+i.detail["path"]+"\n")
                                for j in sub_list:
                                    file.write("sub "+str(j[0])+"\n")
                                    for k in j[1].distance_from_point_to_line(j[1].label[0], )[:20]:
                                        file.write(str(k)+"\n")
                                    if not os.path.isdir(i.detail["path"]+"/"+str(sub_list.index(j))):
                                        linux_command("mkdir "+i.detail["path"]+"/"+str(sub_list.index(j)))
                                        j[1].write(i.detail["path"]+"/"+str(sub_list.index(j))+"/CONTCAR")
                                    i.bundle.append(i.detail["path"]+"/"+str(sub_list.index(j)))
                            i.state="done"
                                
                        else:
                            p.read(i.detail["structure"]+"/CONTCAR")
                            sub_list=p.sub(i.detail["origin_site"],i.detail["new_site"])
                            if not os.path.isdir(i.detail["path"]):
                                linux_command("mkdir "+i.detail["path"])
                            if len(sub_list)==1:
                                sub_list[0][1].write(i.detail["path"]+"/CONTCAR")
                            else:
                                file.write("\n"+i.detail["path"]+"\n")
                                for j in sub_list:
                                    file.write("sub "+str(j[0])+"\n")
                                    for k in j[1].distance(j[1].label[0],True)[:20]:
                                        file.write(str(k)+"\n")
                                    if not os.path.isdir(i.detail["path"]+"/"+j[0][0]+"_"+str(j[0][1])):
                                        linux_command("mkdir "+i.detail["path"]+"/"+j[0][0]+"_"+str(j[0][1]))
                                        j[1].write(i.detail["path"]+"/"+j[0][0]+"_"+str(j[0][1])+"/CONTCAR")
                                    i.bundle.append(i.detail["path"]+"/"+j[0][0]+"_"+str(j[0][1]))
                            i.state="done"

                elif i.detail["type"]=="cc":
                    if os.path.isdir(i.detail["path"]):
                        for j in os.listdir(i.detail["path"]):
                            i.bundle.append(i.detail["path"]+"/"+j)
                        i.state="done"
                    if (self.find_job_by_path(i.detail["structure"]).state=="done" and 
                        self.find_job_by_path(i.detail["vib1"]).state=="done" and
                        self.find_job_by_path(i.detail["vib2"]).state=="done"):
                        p=poscar()
                        v1=poscar()
                        v2=poscar()
                        p.read(i.detail["structure"]+"/CONTCAR")
                        v1.read(i.detail["vib1"]+"/CONTCAR")
                        v2.read(i.detail["vib2"]+"/CONTCAR")
                        if not os.path.isdir(i.detail["path"]):
                            linux_command("mkdir "+i.detail["path"])
                        for j in p.cc(v1-v2,i.detail["ini"],i.detail["fin"],i.detail["step"]):
                            if not os.path.isdir(i.detail["path"]+"/cc_"+str(j[0])):
                                linux_command("mkdir "+i.detail["path"]+"/cc_"+str(j[0]))
                                j[1].write(i.detail["path"]+"/cc_"+str(j[0])+"/CONTCAR")
                            i.bundle.append(i.detail["path"]+"/cc_"+str(j[0]))
                        i.state="done"
                
                elif i.detail["type"]=="cc2":
                    if os.path.isdir(i.detail["path"]):
                        for j in os.listdir(i.detail["path"]):
                            i.bundle.append(i.detail["path"]+"/"+j)
                        i.state="done"
                    if (self.find_job_by_path(i.detail["structure"]).state=="done" and 
                        self.find_job_by_path(i.detail["vib1"]).state=="done" and
                        self.find_job_by_path(i.detail["vib2"]).state=="done"):
                        p=poscar()
                        v1=poscar()
                        v2=poscar()
                        p.read(i.detail["structure"]+"/CONTCAR")
                        v1.read(i.detail["vib1"]+"/CONTCAR")
                        v2.read(i.detail["vib2"]+"/CONTCAR")
                        if not os.path.isdir(i.detail["path"]):
                            linux_command("mkdir "+i.detail["path"])
                        for j in p.cc(v1-v2,i.detail["ini"],i.detail["fin"],i.detail["step"]):
                            for k in j[1].cc(p.cc2(),0,0.5,6):
                                if not os.path.isdir(i.detail["path"]+"/cc_"+str(j[0])+"_"+str(k[0])):
                                    linux_command("mkdir "+i.detail["path"]+"/cc_"+str(j[0])+"_"+str(k[0]))
                                    k[1].write(i.detail["path"]+"/cc_"+str(j[0])+"_"+str(k[0])+"/CONTCAR")
                                i.bundle.append(i.detail["path"]+"/cc_"+str(j[0])+"_"+str(k[0]))
                        i.state="done"
                
                elif i.detail["type"]=="distort":
                    if os.path.isdir(i.detail["path"]):
                        i.state="done"
                    elif i.detail["structure"] in self.job_path_list and self.find_job_by_path(i.detail["structure"]).state=="done":
                        p=poscar()
                        if self.find_job_by_path(i.detail["structure"]).bundle:
                            pass
                        else:
                            p.read(i.detail["structure"]+"/CONTCAR")
                            linux_command("mkdir "+i.detail["path"])
                            p.distort(i.detail["type"]).write(i.detail["path"]+"/CONTCAR")
                        i.state="done"
                
                elif i.detail["type"]=="distort_sphere":
                    if os.path.isdir(i.detail["path"]):
                        i.state="done"
                    elif (i.detail["structure"] in self.job_path_list
                          and self.find_job_by_path(i.detail["structure"]).state=="done"
                          and (not "/" in i.detail["distort_list"]
                               or ("/" in i.detail["distort_list"] 
                                   and self.find_job_by_path(i.detail["distort_list"]).state=="done"
                                   and self.find_job_by_path(origin_name(i.detail["distort_list"])).state=="done"))):
                        
                        if i.detail["distort_list"]=="0":
                            i.state="fail"
                        else:
                            if "/" in i.detail["distort_list"]:
                                distort_poscar=poscar()
                                origin_poscar=poscar()
                                origin_poscar.read(origin_name(i.detail["distort_list"])+"/CONTCAR")
                                distort_poscar.read(i.detail["distort_list"]+"/CONTCAR")
                                if "half" in i.detail["path"]:
                                    distortion=(distort_poscar-origin_poscar)*0.5
                                elif "double" in i.detail["path"]:
                                    distortion=(distort_poscar-origin_poscar)*2
                                else:
                                    distortion=distort_poscar-origin_poscar
                                distort_list=[]
                                current_element=False
                                neis=0
                                for j in distort_poscar.distance(distort_poscar.label[0],False):
                                    if not current_element==j[0][0]:
                                        neis+=1
                                        if neis>3:
                                            break
                                        else:
                                            current_element=j[0][0]
                                    abs_origin=(origin_poscar.position[distort_poscar.label.index(j[0])]
                                            -origin_poscar.position[0]).dot(origin_poscar.axis)
                                    abs_distortion=distortion.position[distort_poscar.label.index(j[0])].dot(distortion.axis)
                                    if np.linalg.norm(abs_origin):
                                        distort_list.append(round(abs_origin.dot(abs_distortion)/np.linalg.norm(abs_origin),3))
                            else:
                                distort_list=[float(j) for j in i.detail["distort_list"].split("_")]
                            high_sym_poscar=poscar()
                            high_sym_poscar.read(i.detail["structure"]+"/CONTCAR")
                            high_sym_distance=high_sym_poscar.distance(high_sym_poscar.label[0],False)
                            if (distort_list[5]-distort_list[0])<1e-5:
                                distort_mode="oh"
                            elif ((distort_list[3]-distort_list[0])<1e-3
                                and (distort_list[5]-distort_list[4])<1e-3):
                                distort_mode="ez"
                            elif ((distort_list[1]-distort_list[0])<1e-3
                                  and (distort_list[5]-distort_list[2])<1e-3):
                                distort_mode="cz"
                            else:
                                distort_mode="nosym"
                            
                            avg=0
                            z_list=[]
                            for j in range(7)[1:]:
                                avg+=high_sym_distance[j][1]
                                pos=(high_sym_poscar.position[high_sym_poscar.label.index(high_sym_distance[j][0])]
                                     -high_sym_poscar.position[0])
                                if abs(pos[0])<1e-3 and abs(pos[1])<1e-3:
                                    z_list.append(j)
                            avg/=6
                            pre_distort_list=[]
                            if distort_mode=="ez":
                                for j in range(7)[1:]:
                                    pre_distort_list.append(avg-high_sym_distance[j][1]+(1e-5 if j in z_list else -1e-5))
                            elif distort_mode=="cz":
                                for j in range(7)[1:]:
                                    pre_distort_list.append(avg-high_sym_distance[j][1]-(1e-5 if j in z_list else -1e-5))
                            linux_command("mkdir "+i.detail["path"])
                            high_sym_poscar.distort_sphere(pre_distort_list).distort_sphere(distort_list).write(i.detail["path"]+"/CONTCAR")
                            i.state="done"
                
                elif i.detail["type"]=="expand":
                    if os.path.isdir(i.detail["path"]):
                        i.state="done"
                    elif self.find_job_by_path(i.detail["structure"]).state=="done":
                        p=poscar()
                        if self.find_job_by_path(i.detail["structure"]).bundle:
                            pass
                        else:
                            p.read(i.detail["structure"]+"/CONTCAR")
                            linux_command("mkdir "+i.detail["path"])
                            p.expand(i.detail["tm"])
                            p.write(i.detail["path"]+"/CONTCAR")
                        i.state="done"
                
                elif "slab" in i.detail["type"]:
                    if os.path.isdir(i.detail["path"]):
                        i.state="done"
                    elif self.find_job_by_path(i.detail["structure"]).state=="done":
                        p=poscar()
                        if self.find_job_by_path(i.detail["structure"]).bundle:
                            pass
                        else:
                            p.read(i.detail["structure"]+"/CONTCAR")
                            linux_command("mkdir "+i.detail["path"])
                            if i.detail["type"]=="slab":
                                p.slab(i.detail["vac_length"])
                            elif i.detail["type"]=="slab2":
                                p.slab2(i.detail["vac_length"])
                            p.write(i.detail["path"]+"/CONTCAR")
                        i.state="done"
                
                elif "primcell" in i.detail["type"]:
                    if os.path.isdir(i.detail["path"]):
                        i.state="done"
                    elif self.find_job_by_path(i.detail["structure"]).state=="done":
                        linux_command("mkdir "+i.detail["path"])
                        linux_command("cp "+i.detail["structure"]+"/CONTCAR "+i.detail["path"]+"/POSCAR")
                        cwd=os.getcwd()
                        os.chdir(i.detail["path"])
                        linux_command("vaspkit -task 303")
                        os.chdir(cwd)
                        linux_command("mv "+i.detail["path"]+"/PRIMCELL.vasp "+i.detail["path"]+"/CONTCAR")
                        i.state="done"
                
                elif i.detail["type"]=="plot_pot":
                    if False and os.path.isfile(i.detail["path"]+".png"):
                        i.state="done"
                    elif self.find_job_by_path(i.detail["pot"]).state=="done":
                        procar=read_procar(i.detail["pot"])[0]
                        for j in procar:
                            if j["occ"]<0.3:
                                vbm=procar[procar.index(j)-1]["energy"]
                                break
                        #file.write("\n"+i.detail["path"]+" "+str(vbm)+" ")
                        cwd=os.getcwd()
                        os.chdir(i.detail["pot"])
                        linux_command("(echo 426;echo 3)|vaspkit")
                        os.chdir(cwd)
                        
                        with open(i.detail["pot"]+"/PLANAR_AVERAGE.dat","r") as f:
                            lines=f.readlines()
                        match_result=re.search("#(.+), (.+)",lines[0])
                        if not match_result:
                            match_result=re.search("#(.+) (.+)",lines[0])
                        xlabel,ylabel=match_result.groups()
                        x=[xlabel]
                        y=[ylabel]
                        for j in lines[1:]:
                            count=0
                            for k in j.replace("\n","").split(" "):
                                if k:
                                    count+=1
                                    if count%2==1:
                                        x.append(float(k))
                                    else:
                                        y.append(float(k))
                        
                        count={}
                        for j in [round(j,3) for j in y[1:]]:
                            if j in count:
                                count[j]+=1
                            else:
                                count[j]=1
                        
                        count=1
                        ysum=y[1]
                        temp_count=[]
                        temp_ysum=[]
                        file.write("avg "+i.detail["pot"]+"\n")
                        for j in range(len(x))[2:-1]:
                            count+=1
                            ysum+=y[j]
                            if (y[j]-y[j-1])/(x[j]-x[j-1])*(y[j+1]-y[j])/(x[j+1]-x[j])<0:
                                temp_count.append(count)
                                temp_ysum.append(ysum)
                                #file.write(str(count)+" "+str(ysum)+" ")
                                count=0
                                ysum=0
                        
                        for j in temp_count:
                            file.write(str(j)+" ")
                        file.write("\n")
                        for j in temp_ysum:
                            file.write(str(j)+" ")
                        file.write("\n")
                        for j,k in zip(temp_count,temp_ysum):
                            file.write(str(k/j)+" ")
                        file.write("\n")
                        
                        count=0
                        ysum=0
                        file.write("sum "+i.detail["pot"]+" ")
                        for j in range(len(x))[1:]:
                            count+=1
                            ysum+=y[j]
                            file.write(str(ysum/count)+" ")
                        file.write("\n")
                        
                        good_range=3
                        good_list=[]
                        for j in range(len(x))[good_range+1:len(x)-good_range]:
                            if abs((y[j]-y[j-good_range])/(x[j]-x[j-good_range])-
                                   (y[j]-y[j+good_range])/(x[j]-x[j+good_range]))<0.01:
                                good_list.append([x[j],y[j]])
                        if len(good_list)>5:
                            file.write("pot "+
                                       i.detail["pot"]+" "+
                                       str(sum(list(zip(*good_list))[1])/len(good_list))+" "+
                                       str(good_list[int(len(good_list)/2)][1])+"\n")
                        else:
                            # file.write("pot "+
                            #            i.detail["pot"]+" "+
                            #            str(y[1])+"\n")
                            pass
                        
                        
                        plot_pot(x,y,i.detail["path"])
                        i.state="done"
                        
                elif i.detail["type"]=="plot_cc":
                    if not i.bundle:
                        # vars in i.detail
                        for j in ["grd","ex","ex_at_grd","grd_at_ex"]:
                            i.bundle.append(self.find_job_by_path(i.detail[j]))
                    if set([j.state for j in i.bundle])=={"done"}:
                        file.write("\n"+i.detail["path"]+"\n")
                        
                        energy_list=[energy(j.detail["path"]) if j else j for j in i.bundle]
                        if not os.path.isfile(i.detail["path"]+".png"):
                            plot_cc(energy_list,i.detail["path"])
                        i.state="done"
                elif i.detail["type"]=="plot_cc2" and not delay:
                    energy_dict1={}
                    system=i.detail["path"].split("/")[0]
                    p=poscar()
                    for j in self.job_path_list:
                        if j in [system+"/"+k 
                                 for k in ["ex_at_grd_pbe0","ex_at_grd_pbe0_soc","ex_cb_at_grd_pbe0_soc","ex_vb_at_grd_pbe0_soc",
                                           "ex_ez_pbe0_soc_temp1",
                                           "grd_at_ex_cz_pbe0_soc","grd_at_ex_ez_pbe0","grd_at_ex_ez_pbe0_soc","grd_at_ex_cb_pbe0_soc","grd_at_ex_vb_pbe0_soc",
                                           "relax_ex_cb_pbe0_soc","relax_ex_cz_pbe0_soc","relax_ex_ez_pbe0","relax_ex_ez_pbe0_soc","relax_ex_vb_pbe0_soc",
                                           "relax_grd_pbe0","relax_grd_pbe0_soc",
                                           "relax_grdm1_pbe0_soc",
                                           "relax_grdp1_pbe0_soc"]] and self.find_job_by_path(j).state=="done":
                            p.read(j+"/CONTCAR")
                            energy_dict1[j.split("/")[-1]]=(p.cc_xtick(i.detail["type"],5 if "vcl" in system else 6),energy(j))
                        elif j in [system+"/grdm1_pbe0_soc",
                                   system+"/relax_ex_cz_half_pbe0_soc",
                                   system+"/relax_ex_ez_half_pbe0_soc",
                                   system.replace("1","0").split("_")[0]+"/unitcell_pbe0_soc"] and self.find_job_by_path(j).state=="done":
                                       
                            procar=read_procar(j)[0]
                            sub_atom=j.split("/")[0].split("_")[-1]
                            if j.split("/")[-1]=="grdm1_pbe0_soc":
                                orbit_dict={"vb":False,
                                           "s":False,
                                           "p":False,
                                           "cb":False}
                                for k in procar:
                                    if k["occ"]<=0.5:
                                        if (sub_atom,"s") in k["spd"] and k["spd"][(sub_atom,"s")]>0.05:
                                            if not orbit_dict["s"]:
                                                orbit_dict["s"]=k["energy"]
                                        elif (sub_atom,"p") in k["spd"] and k["spd"][(sub_atom,"p")]>0.1:
                                            if not orbit_dict["p"]:
                                                orbit_dict["p"]=k["energy"]
                                        elif (not k["spd"] or not main_comp(k["spd"])[0] in ["F","Cl","Br"]):
                                            if not orbit_dict["cb"]:
                                                orbit_dict["cb"]=k["energy"]
                                        elif not orbit_dict["p"] and main_comp(k["spd"])[0] in ["F","Cl","Br"]:
                                            if not orbit_dict["vb"]:
                                                orbit_dict["vb"]=k["energy"]
                                    if orbit_dict["p"] and orbit_dict["cb"]:
                                        
                                        break
                                energy_dict1[j.split("/")[-1]]=orbit_dict
                                #energy_dict1[j.split("/")[-1]]=dict(list(orbit_dict.items())[:])
                            elif j.split("/")[-1]=="unitcell_pbe0_soc":
                                energy_dict1[j.split("/")[-1]]={}
                                for k in procar:
                                    if k["occ"]<0.5:
                                        energy_dict1[j.split("/")[-1]]["cbm"]=k["energy"]
                                        energy_dict1[j.split("/")[-1]]["vbm"]=procar[procar.index(k)-1]["energy"]
                                        break
                            else:
                                p.read(j+"/CONTCAR")
                                energy_dict1[j.split("/")[-1]]=(p.cc_xtick(i.detail["type"],5 if "vcl" in system else 6),{})
                                for k in procar:
                                    if k["occ"]<0.4:
                                        break
                                    elif k["occ"]<0.6:
                                        if (sub_atom,"s") in k["spd"] and not "s" in energy_dict1[j.split("/")[-1]][1]:
                                            energy_dict1[j.split("/")[-1]][1]["s"]=k["energy"]
                                        elif (sub_atom,"p") in k["spd"] and not "p" in energy_dict1[j.split("/")[-1]][1]:
                                            energy_dict1[j.split("/")[-1]][1]["p"]=k["energy"]
                    
                    energy_dict1["ex_vb_at_grd_m1_pbe0_soc"]=False
                    energy_dict1["ex_vc_at_grd_m1_pbe0_soc"]=False
                    energy_dict1["ex_at_grd_m1_pbe0_soc"]=False
                    energy_dict1["ex_cb_at_grd_m1_pbe0_soc"]=False
                    for j in energy_dict1:
                        if j in ["ex_at_grd_pbe0",
                                 "grd_at_ex_ez_pbe0",
                                 "relax_ex_ez_pbe0"]:
                            energy_dict1[j]=(energy_dict1[j][0],
                                             energy_dict1[j][1]-energy_dict1["relax_grd_pbe0"][1]+energy_dict1["relax_grd_pbe0_soc"][1])
                        elif j=="relax_grdm1_pbe0_soc":
                            energy_dict1[j]=(energy_dict1[j][0],
                                             energy_dict1[j][1]+energy_dict1["unitcell_pbe0_soc"]["cbm"])
                        elif j=="relax_grdp1_pbe0_soc":
                            energy_dict1[j]=(energy_dict1[j][0],
                                             energy_dict1[j][1]-energy_dict1["unitcell_pbe0_soc"]["vbm"])
                        elif j in ["relax_ex_cz_half_pbe0_soc","relax_ex_ez_half_pbe0_soc"]:
                            energy_dict1[j]=(energy_dict1[j][0],
                                             energy_dict1[j][1]["p"]-energy_dict1[j][1]["s"]+energy_dict1["relax_grd_pbe0_soc"][1])
                        elif j=="grdm1_pbe0_soc":
                            if energy_dict1[j]["vb"]:
                                if energy_dict1[j]["p"]:
                                    energy_dict1["ex_vb_at_grd_m1_pbe0_soc"]=(energy_dict1["relax_grd_pbe0_soc"][0],
                                                                              energy_dict1[j]["p"]-energy_dict1[j]["vb"]+energy_dict1["relax_grd_pbe0_soc"][1])
                                if energy_dict1[j]["cb"]:
                                    energy_dict1["ex_vc_at_grd_m1_pbe0_soc"]=(energy_dict1["relax_grd_pbe0_soc"][0],
                                                                              energy_dict1[j]["cb"]-energy_dict1[j]["vb"]+energy_dict1["relax_grd_pbe0_soc"][1])
                            if energy_dict1[j]["s"]:
                                if energy_dict1[j]["p"]:
                                    energy_dict1["ex_at_grd_m1_pbe0_soc"]=(energy_dict1["relax_grd_pbe0_soc"][0],
                                                                           energy_dict1[j]["p"]-energy_dict1[j]["s"]+energy_dict1["relax_grd_pbe0_soc"][1])
                                if energy_dict1[j]["cb"]:
                                    energy_dict1["ex_cb_at_grd_m1_pbe0_soc"]=(energy_dict1["relax_grd_pbe0_soc"][0],
                                                                              energy_dict1[j]["cb"]-energy_dict1[j]["s"]+energy_dict1["relax_grd_pbe0_soc"][1])
                    
                    # if "ex_at_grd_pbe0_soc" in energy_dict1:
                    #     if "ex_at_grd_m1_pbe0_soc" in energy_dict1:
                    #         del energy_dict1["ex_at_grd_m1_pbe0_soc"]
                    #     if "relax_ex_ez_pbe0_soc" in energy_dict1:
                    #         for j in ["ex_at_grd_pbe0",
                    #                   "ex_ez_pbe0_soc_temp1",
                    #                   "grd_at_ex_ez_pbe0",
                    #                   "relax_ex_ez_pbe0",
                    #                   "relax_grd_pbe0"]:
                    #             if j in energy_dict1:
                    #                 del energy_dict1[j]
                    
                    if "relax_grd_pbe0" in energy_dict1 and "relax_grd_pbe0_soc" in energy_dict1:
                        energy_dict1["relax_grd_pbe0"]=(energy_dict1["relax_grd_pbe0"][0],
                                                    energy_dict1["relax_grd_pbe0_soc"][1])
                    
                    energy_dict2={}
                    distance_dict={}
                    for j in energy_dict1:
                        if type(energy_dict1[j])==tuple:
                            energy_dict2[replace_dict[j]]=energy_dict1[j][1]
                            distance_dict[replace_dict[j]]=energy_dict1[j][0]
                    energy_dict2["system"]=system
                    distance_dict["system"]=system
                    #file.write(str(energy_dict1)+"ccc\n")
                    if i.detail["type"]=="dist":
                        file.write(str(energy_dict2)+"cce\n")
                    #file.write(str(distance_dict)+"ccd\n")
                    file.write(system+"ccc"+str(energy_dict1)+"\n")
                    if not os.path.isfile(i.detail["path"]+".png"):
                        plot_cc2(energy_dict1,i.detail["path"])
                    i.state="done"
                        
                elif i.detail["type"]=="plot_split":
                    if not i.bundle:
                        # vars in i.detail
                        for j in ["exm1","exm1_soc","grd_soc"]:
                            i.bundle.append(self.find_job_by_path(i.detail[j]))
                    if set([j.state for j in i.bundle])=={"done"}:
                        spd=tuple(i.detail["spd"].split("_"))
                        split_list=[]
                        for j in i.bundle:
                            # procar_list结构决定
                            current_split=[]
                            for k in read_procar(j.detail["path"]):
                                for l in k:
                                    # 应对简并轨道,必要时降低标准
                                    if spd in l["spd"] and l["spd"][spd]>0.09 and l["occ"]<0.3:
                                        current_split.append(l)
                            dege={
                                "s":2,
                                "p":6,
                                "d":10,
                                "f":14,
                                }
                            # split_list.append(sorted([m["energy"] for m in sorted(current_split,key=lambda x:x["spd"][spd],reverse=True)[:dege["p"]]]))
                            current_split.sort(key=lambda x:x["energy"])
                            split_list.append([m["energy"] for m in current_split[:dege["p"]]])
                        file.write("\n")
                        for j,l in zip(split_list,["1cfs","2cfs&sos","3sos"]):
                            file.write(l+" [")
                            for k in j:
                                file.write(str(k)+",")
                            file.write("\""+i.detail["path"]+"\"],\n")
                        if not os.path.isfile(i.detail["path"]+".png"):
                            plot_split(split_list,i.detail["path"])
                        i.state="done"
                
                elif i.detail["type"]=="plot_split_cc":
                    if self.find_job_by_path(i.detail["cc_path"]).state=="done":
                        file.write("\n"+i.detail["path"]+"\n")
                        spd=tuple(i.detail["spd"].split("_"))
                        split_list=[]
                        for j in self.find_job_by_path(i.detail["cc_path"]).bundle:
                            # procar_list结构决定
                            current_split=[]
                            for k in read_procar(j):
                                for l in k:
                                    # 应对简并轨道,必要时降低标准
                                    if spd in l["spd"] and l["spd"][spd]>0.09 and l["occ"]<0.3:
                                        current_split.append(l)
                            dege={
                                "s":2,
                                "p":6,
                                "d":10,
                                "f":14,
                                }
                            current_split.sort(key=lambda x:x["energy"])
                            split_list.append([[m["energy"] for m in current_split[:dege["p"]]],float(j.split("/")[-1].replace("cc_",""))])
                        split_list.sort(key=lambda x:x[1])
                        for j in split_list:
                            file.write(str(j)+"\n")
                        if not os.path.isfile(i.detail["path"]+".png"):
                            plot_split_cc(split_list,i.detail["path"])
                        i.state="done"
                
                elif i.detail["type"]=="plot_split_cc2":
                    if (self.find_job_by_path(i.detail["cc1_path"]).state=="done" 
                        and self.find_job_by_path(i.detail["cc2_path"]).state=="done"):
                        spd=tuple(i.detail["spd"].split("_"))
                        split_list=[]
                        for j in self.find_job_by_path(i.detail["cc1_path"]).bundle:
                            # procar_list结构决定
                            current_split=[]
                            for k in read_procar(j):
                                for l in k:
                                    # 应对简并轨道,必要时降低标准
                                    if spd in l["spd"] and l["spd"][spd]>0.09 and l["occ"]<0.3:
                                        current_split.append(l)
                            dege={
                                "s":2,
                                "p":6,
                                "d":10,
                                "f":14,
                                }
                            current_split.sort(key=lambda x:x["energy"])
                            split_list.append([[m["energy"] for m in current_split[:dege["p"]]],-float(j.split("/")[-1].replace("cc_",""))])
                        
                        for j in self.find_job_by_path(i.detail["cc2_path"]).bundle:
                            # procar_list结构决定
                            current_split=[]
                            for k in read_procar(j):
                                for l in k:
                                    # 应对简并轨道,必要时降低标准
                                    if spd in l["spd"] and l["spd"][spd]>0.09 and l["occ"]<0.3:
                                        current_split.append(l)
                            dege={
                                "s":2,
                                "p":6,
                                "d":10,
                                "f":14,
                                }
                            current_split.sort(key=lambda x:x["energy"])
                            split_list.append([[m["energy"] for m in current_split[:dege["p"]]],float(j.split("/")[-1].replace("cc_",""))])
                        split_list.sort(key=lambda x:x[1])
                        if not os.path.isfile(i.detail["path"]+".png"):
                            plot_split_cc(split_list,i.detail["path"])
                        i.state="done"
                
                elif i.detail["type"]=="plot_ctl":
                    if not i.bundle:
                        ctl_pattern=re.compile(i.detail["pattern"].replace("xxx","(\w)(\d)")+"$")
                        for j in self.job_list:
                            if "path" in j.detail:
                                match_result=ctl_pattern.search(j.detail["path"])
                                if match_result:
                                    i.bundle.append(j)
                        i.bundle.append(self.find_job_by_path(i.detail["pattern"].replace("xxx","")))
                        i.bundle.append(self.find_job_by_path(i.detail["unitcell"]))
                    #if set([j.state for j in i.bundle])=={"done"}:
                    if i.bundle[-1].state=="done":
                        # bandgap
                        job_unitcell=self.find_job_by_path(i.detail["unitcell"])
                        procar=read_procar(job_unitcell.detail["path"])[0]
                        for j in procar:
                            if j["occ"]<0.3:
                                cbm=j["energy"]
                                vbm=procar[procar.index(j)-1]["energy"]
                                break
                        energy_unitcell=energy(job_unitcell.detail["path"])
                        # 是超胞
                        if "supercell" in self.find_job_by_path(i.detail["path"].split("/")[0]+"/sub").detail["structure"]:
                            expand_job=self.find_job_by_path("0"+i.detail["path"].split("/")[0].split("_")[0][1:]+"/supercell")
                            energy_unitcell*=np.linalg.det(expand_job.detail["tm"])
                        energy_list=[[-j.detail["nele"],energy(j.detail["path"])-energy_unitcell] for j in i.bundle[:-1] if j.state=="done"]
                        
                        energy_list.sort(key=lambda x:x[0])
                        file.write("\n["+str(vbm)+","+str(cbm)+",")
                        for j in range(len(energy_list)-1):
                            file.write(str(energy_list[j][1]-energy_list[j+1][1])+",")
                        file.write("\""+i.detail["path"]+"\"],\n")
                        for j in energy_list:
                            file.write(str(j[0])+" "+str(j[1])+"\n")
                        if not os.path.isfile(i.detail["path"]+".png"):
                            plot_ctl([vbm,cbm],energy_list,i.detail["path"])
                        i.state="done"
                
                elif i.detail["type"]=="plot_ctl3":
                    if not i.bundle:
                        for j in [i.detail["p1"],
                                  i.detail["p1"].replace("p1",""),
                                  i.detail["m1"],
                                  i.detail["m1"].replace("m1",""),
                                  i.detail["unitcell"]]:
                            i.bundle.append(self.find_job_by_path(j))
                    
                    if set([j.state for j in i.bundle])in [{"done","wrong"},{"done"}]:
                        # bandgap
                        procar=read_procar(i.bundle[-1].detail["path"])[0]
                        for j in procar:
                            if j["occ"]<0.3:
                                cbm=j["energy"]
                                vbm=procar[procar.index(j)-1]["energy"]
                                break
                        
                        energy_list=[[-i.bundle[0].detail["nele"],energy(i.bundle[0].detail["path"])-energy(i.bundle[1].detail["path"])],
                                     [-i.bundle[1].detail["nele"],0],
                                     [-i.bundle[2].detail["nele"],energy(i.bundle[2].detail["path"])-energy(i.bundle[3].detail["path"])]]
                        
                        energy_list.sort(key=lambda x:x[0])
                        file.write("\n["+str(vbm)+","+str(cbm)+",")
                        for j in range(len(energy_list)-1):
                            if i.bundle[j].state=="wrong":
                                file.write(str(1000+energy_list[j][1]-energy_list[j+1][1])+",")
                            else:
                                file.write(str(energy_list[j][1]-energy_list[j+1][1])+",")
                        file.write("\""+i.detail["path"]+"\"],\n")
                        for j in energy_list:
                            file.write(str(j[0])+" "+str(j[1])+"\n")
                        if not os.path.isfile(i.detail["path"]+".png"):
                            plot_ctl([vbm,cbm],energy_list,i.detail["path"])
                        i.state="done"
                
                elif i.detail["type"]=="ctl4" and not delay:

                    if self.find_job_by_path(i.detail["system"].replace("1","0").split("_")[0]+"/unitcell_"+i.detail["level"]).state=="done":
                        energy_dict1={}
                        energy_dict2={}
                        for j in self.job_path_list:
                            for k in ["grdm1","grdp1","ex","grd"]:
                                if i.detail["system"]+"/relax_"+k in j and i.detail["level"] in j:
                                    if self.find_job_by_path(j).state=="done":
                                        energy_dict1[j.replace(i.detail["system"]+"/relax_","").replace("_"+i.detail["level"],"")]=energy(j)-charge_corr_dict[i.detail["system"].replace("1","0").split("_")[0]]*self.find_job_by_path(j).detail["nele"]**2
                                    elif self.find_job_by_path(j).state=="wrong":
                                        energy_dict1[j.replace(i.detail["system"]+"/relax_","").replace("_"+i.detail["level"],"")]=energy(j)+1e5-charge_corr_dict[i.detail["system"].replace("1","0").split("_")[0]]*self.find_job_by_path(j).detail["nele"]**2
                                    break
                        
                        if "grd" in energy_dict1:
                            for j in energy_dict1:
                                if "grdp1" in j:
                                    energy_dict2["0/1"+j.replace("grdp1","")]=energy_dict1[j]-energy_dict1["grd"]
                        
                        if "grdm1" in energy_dict1:
                            for j in energy_dict1:
                                if not "grdp1" in j and not "grdm1" in j:
                                    energy_dict2["-1/0"+j.replace("grd","").replace("ex","")]=energy_dict1[j]-energy_dict1["grdm1"]
                        
                        procar=read_procar(i.detail["system"].replace("1","0").split("_")[0]+"/unitcell_"+i.detail["level"])[0]
                        for j in procar:
                            if j["occ"]<0.3:
                                energy_dict2["CBM"]=j["energy"]
                                energy_dict2["VBM"]=procar[procar.index(j)-1]["energy"]
                                break
                        
                        energy_dict2["system"]=i.detail["system"]
                        
                        file.write(str(energy_dict2)+"ctl4\n")
                        i.state="done"    
                        
                # update plot_nei要在所有bundle的vasp计算之后
                elif i.detail["type"]=="plot_nei" and not delay:
                    for j in self.job_list:
                        if ("path" in j.detail
                            and j.detail["path"].split("/")[0]==i.detail["path"].split("/")[0]
                            and "relax" in j.detail["path"]
                            and j.state=="done"
                            and not j.bundle):
                            i.bundle.append(j)
                    # 多于一个才可以画子图
                    if len(i.bundle)>0:
                        file.write("\n"+i.detail["path"]+"\n")
                        p=poscar()
                        nei_list=[]
                        match_result=re.search("(\w+)_(\d+)",i.detail["atom"])
                        for j in i.bundle:
                            p.read(j.detail["path"]+"/CONTCAR")
                            current_nei=[]
                            for k in p.env((match_result.group(1),int(match_result.group(2))),2,False):
                                if k[1]<i.detail["range"]:
                                    current_nei.append(k)
                                else:
                                    current_nei.append(j.detail["path"])
                                    break
                            nei_list.append(current_nei)
                        for j in nei_list:
                            file.write(str(j)+",\n")
                        if len(i.bundle)>1:
                            plot_nei(nei_list,i.detail["path"])
                        i.state="done"
                        
                elif i.detail["type"]=="plot_band":
                    if os.path.isfile(i.detail["path"]+".png"):
                        i.state="done"
                    elif self.find_job_by_path(i.detail["structure"]).state=="done":
                        cwd=os.getcwd()
                        os.chdir(i.detail["structure"])
                        linux_command("vaspkit -task 211")
                        band_list=[]
                        label=[[],[]]
                        with open("REFORMATTED_BAND_UP.dat","r") as f:
                            lines=f.readlines()
                            for j in lines[1:]:
                                temp_path=[float(k) for k in j.split(" ") if not k==""]
                                if band_list:
                                    for k in range(len(temp_path)):
                                        band_list[k].append(temp_path[k])
                                else:
                                    for k in temp_path:
                                        band_list.append([k])
                                        
                        with open("KLABELS","r") as f:
                            lines=f.readlines()
                            for j in lines[1:]:
                                temp_label=[float(k) if "." in k else k for k in j.split(" ") if not k==""]
                                if len(temp_label)<2:
                                    break
                                label[0].append(temp_label[1])
                                label[1].append(temp_label[0])
                        os.chdir(cwd)
                        plot_band(band_list,label,i.detail["path"])
                        i.state=="done"

        return None
    
    def run(self,file):
        import os
        import os.path
        import time
        
    
        for i in list(set([i.split("/")[0] for i in self.job_path_list])):
            if not os.path.isdir(i):
                linux_command("mkdir "+i)
            
        # 清理过时任务
        # import re
        # jobid_template=re.compile("(\d+)")
        # for i in os.listdir():
        #     if os.path.isdir(i):
        #         for j in os.listdir(i):
        #             if os.path.isdir(i+"/"+j):
        #                 listdir=os.listdir(i+"/"+j)
        #                 if "CONTCAR" in listdir and "PROCAR" in listdir and "jobid" in listdir and i[0] in "123456789" and "OSZICAR" in listdir:
        #                     with open(i+"/"+j+"/jobid","r") as f:
        #                         jobid=jobid_template.search(f.readlines()[-1]).group(1)
        #                         print(i+"/"+j,end=" ")
        #                     if jobid+".log" in listdir:
        #                         if "relax" in j:
        #                             with open(i+"/"+j+"/"+jobid+".log","r") as f:
        #                                 log_lines=f.readlines()
        #                             if "for stderr output of this job" in log_lines[-2]:
        #                                 log_lines=log_lines[:-6]
        #                             if "reached required accuracy - stopping structural energy minimisation" in log_lines[-2]:
        #                                 p=poscar()
        #                                 p.read(i+"/"+j+"/CONTCAR")
        #                                 print(p.distance(p.label[0],False)[:7],end=" ")
        #                                 print(energy(i+"/"+j),end=" ")
        #                                 procar=read_procar(i+"/"+j)
        #                                 spd1=(i.split("_")[1],"s")
        #                                 for k in procar[-1]:
        #                                     if k["occ"]<0.7:
        #                                         if not spd1 in k["spd"]:
        #                                             print("ws",end=" ")
        #                                         elif k["occ"]>0.3:
        #                                             print("ws",end=" ")
        #                                         break
                                        
        #                                 spd2=spd1=(i.split("_")[1],"p")
        #                                 for k in procar[0][::-1]:
        #                                     if k["occ"]>0.3:
        #                                         if not spd2 in k["spd"]:
        #                                             print("wp",end=" ")
        #                                         elif k["occ"]<0.7:
        #                                             print("wp",end=" ")
        #                                         break
        #                                 print()
        #                             else:
        #                                 print("notreach")
        #                         else:
        #                             p=poscar()
        #                             p.read(i+"/"+j+"/CONTCAR")
        #                             print(p.distance(p.label[0],False)[:7],end=" ")
        #                             print(energy(i+"/"+j),end=" ")
        #                             procar=read_procar(i+"/"+j)
        #                             spd1=(i.split("_")[1],"s")
        #                             for k in procar[-1]:
        #                                 if k["occ"]<0.7:
        #                                     if not spd1 in k["spd"]:
        #                                         print("ws",end=" ")
        #                                     elif k["occ"]>0.3:
        #                                         print("ws",end=" ")
        #                                     break
                                    
        #                             spd2=spd1=(i.split("_")[1],"p")
        #                             for k in procar[0][::-1]:
        #                                 if k["occ"]>0.3:
        #                                     if not spd2 in k["spd"]:
        #                                         print("wp",end=" ")
        #                                     elif k["occ"]<0.7:
        #                                         print("wp",end=" ")
        #                                     break
        #                             print()
        #                     else:
        #                         print("nolog")
                            
                        # if not i+"/"+j in self.job_path_list:
                        #     if not "notinlist" in i+"/"+j:
                        #         if "jobid" in listdir:
                        #             with open(i+"/"+j+"/jobid","r") as f:
                        #                 linux_command("bkill "+jobid_template.search(f.readlines()[-1]).group(1))
                        #         linux_command("mv "+i+"/"+j+" "+i+"/"+j+"_nil")
                        #     file.write(i+"/"+j+"\n")

                        
        self.update(file,True)
        gonna_finish=False
        while True:
            job_left={}
            job_left["wait"]=[]
            job_left["run"]=[]
            job_left["done"]=[]
            job_left["fail"]=[]
            job_left["wrong"]=[]
            job_dict_dict={}
            
            self.update(file,False)
            for i in self.job_list:
                if "path" in i.detail:
                    job_left[i.state].append(i.detail["path"])
                    if i.detail["path"].split("/")[0][0] in ["1","2"] and i.detail["type"]=="vasp":
                        if not i.detail["path"].split("/")[0] in job_dict_dict:
                            job_dict_dict[i.detail["path"].split("/")[0]]={}
                            job_dict_dict[i.detail["path"].split("/")[0]]["system"]=i.detail["path"].split("/")[0]
                        job_dict_dict[i.detail["path"].split("/")[0]][i.detail["path"].split("/")[1]]=i.state
                        
            with open("state","w") as f:
                print(job_dict_dict.values())
                f.write(str(table(list(job_dict_dict.values())))+"\n")
            
            if not job_left["run"]:
                if gonna_finish:
                    break
                else:
                    gonna_finish=True
                    for i in ["run","fail","wrong"]:
                        print(str(len(job_left[i]))+" "+i)
                        for j in job_left[i]:
                            print(j)
                        print()
            else:
                gonna_finish=False
                for i in ["run","fail","wrong"]:
                    print(str(len(job_left[i]))+" "+i)
                    for j in job_left[i]:
                        print(j)
                    print()
            time.sleep(30)
        return None

def table(dict_list):
    key=[]
    table=[]
    for i in dict_list:
        for j in i.keys():
            if not j in key:
                key.append(j)
    for i in dict_list:
        for j in key:
            if not j in i:
                i[j]=None
    table.append(key)
    for i in dict_list:
        table.append([i[j] for j in key])
    
    for i in table:
        for j in i:
            print(j,end=" ")
        print()
    return None

if __name__ == '__main__':
    with open("report","w") as f:
        jobs().run(f)
        
