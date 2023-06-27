# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 21:48:40 2021

@author: dugue
"""

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
    
    def cc_xtick(self,mode):
        distance_list=[]
        for i in self.distance(self.label[0],False)[1:]:
            if len(distance_list)<6:
                distance_list.append(i[1])
        
        avg=sum(distance_list)/6
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
    
    def slab(self,vac_length):
        import numpy as np
        # z方向扩展2倍
        self.number=[2*i for i in self.number]
        c_length=np.linalg.norm(self.axis[-1])
        ratio=(c_length*2+vac_length)/c_length
        for i in range(len(self.axis[-1])):
            self.axis[-1][i]*=ratio
        new_pos=[]
        for i in range(len(self.position)):
            new_pos.append(np.array([self.position[i][j]/ratio if j==2 else self.position[i][j] for j in range(3)]))
            new_pos.append(np.array([(self.position[i][j]+1)/ratio if j==2 else self.position[i][j] for j in range(3)]))
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

def zsort(poscar_path):
    p=poscar()
    p.read(poscar_path)
    zlist=[]
    for i,j in zip(p.label,p.position):
        if not i[0]=="O":
            zlist.append([len(zlist),j[2]])
    zlist.sort(key=lambda x:x[1])
    print(list(zip(*zlist))[0])
    return None

zsort("ppPOSCAR")

        