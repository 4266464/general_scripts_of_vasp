#!/home/phys/qif/anaconda3/bin/python

import sys
import os

def linux_command(command):
    import os
    with os.popen(command, "r") as p:
        command_return=p.read()
    return command_return

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
            
class poscar:
    def __init__(self):
        pass
    
    def read(self,path):
        import re
        import numpy as np
        scale_pattern=re.compile("([\d\.]+)")
        axis_pattern=re.compile('([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)')
        position_pattern=re.compile('([\d\-\.e]+)\s+([\d\-\.e]+)\s+([\d\-\.e]+)')
        with open(path,"r") as f:
            lines=f.readlines()
        self.comment=lines[0]
        self.scale=float(scale_pattern.search(lines[1]).group(1))
        self.axis=np.array([[float(j) for j in axis_pattern.search(lines[i]).groups()] for i in range(5)[2:]])
        self.element=match_times("\s","\w",lines[5])
        self.number=match_times("\s","\d",lines[6])
        self.label=[(i,k+1) for i,j in zip(self.element,self.number) for k in range(j)]
        self.mode=lines[7]
        self.position=[]
        for i in range(len(lines))[8:]:
            match_result=position_pattern.search(lines[i])
            if match_result:
                self.position.append(np.array([float(j) for j in match_result.groups()]))
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
    
    def distance(self,site,full):
        import numpy as np
        site_position=self.position[self.label.index(site)]
        distance_list=[]
        for i,j in zip(self.label,self.position):
            # 距离
            if full:
                distance_list.append([i,np.array([i if abs(i)<0.5 else (i+1 if i<0 else i-1) for i in j-site_position]).dot(self.axis)])
            else:
                distance_list.append([i,np.linalg.norm(np.array([i if i<0.5 else 1-i for i in abs(site_position-j)]).dot(self.axis))])
        distance_list.sort(key=lambda x:x[0])
        distance_list.sort(key=lambda x:np.linalg.norm(x[1]))
        return distance_list

def read_procar(path):
    import re
    band_pattern=re.compile("band\s+(\d+) # energy\s+([\d\.\-]+) # occ.\s+([\d\.]+)")
    p=poscar()
    p.read(path+"/POSCAR")
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

if __name__ == '__main__':
    if "-e" in sys.argv:
        current_dir=os.getcwd()
        for i,j,k in os.walk("."):
            for l in k:
                if l=="PROCAR":
                    print(i)
                    os.chdir(i)
                    procar=read_procar(".")
                    num=0
                    for m in range(len(procar[0])):
                        if procar[0][m]["occ"]<0.1:
                            if not num:
                                print(procar[0][m])
                                num+=1
                            else:
                                for n in [("Sb","p"),("In","s")]
                            
                            num+=1
                            
                            
                                num=m
                                break
                    for m in range(len(procar[0])):
                        for n in procar:
                            if abs(m-num)<5:
                                print(n[m])
                    os.chdir(current_dir)
    else:
        if len(sys.argv)==1:
            path="."
            num=-1
        elif len(sys.argv)==2:
            if sys.argv[1].isdigit():
                path="."
                num=int(sys.argv[1])
            else:
                path=sys.argv[1]
                num=-1
        else:
            path=sys.argv[1]
            num=int(sys.argv[2])
        
        procar=read_procar(path)
        if num<0:
            for i in range(len(procar[0])):
                if procar[0][i]["occ"]<0.5:
                    num=i
                    break
        for i in range(len(procar[0])):
            for j in procar:
                if abs(i-num)<10:
                    print(j[i])
            
