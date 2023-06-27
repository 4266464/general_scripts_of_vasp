import scipy.special
import math
import numpy as np
from numpy import linalg as la

def nei(radius,x,y,z):
    list_nei=[]
    for i in range(-radius,radius+1):
        for j in range(-radius,radius+1):
            for k in range(-radius,radius+1):
                vec=i*x+j*y+k*z
                dis=round(la.norm(vec),8)
                if not dis > radius*min(la.norm(x),la.norm(y),la.norm(z)) and dis:
                    list_nei.append(vec)
    return list_nei

def ewald_sum(gamma,radius,v,vol,rv,dt,sv):
    res=0
    list_nei_real=nei(radius,v[0],v[1],v[2])
    if np.any(sv):
        list_nei_real.append(np.mat([0,0,0]))
        for i in range(len(list_nei_real)):
            list_nei_real[i]=list_nei_real[i]-sv
    list_nei_rec=nei(radius,rv[0],rv[1],rv[2])
    for i in list_nei_real:
        med=math.sqrt(i*la.inv(dt)*np.transpose(i))
        res += scipy.special.erfc(gamma*med)/med/math.sqrt(la.det(dt))
    for i in list_nei_rec:
        med=float(i*dt*np.transpose(i))
        res += 4*math.pi*np.exp(-med/4/gamma/gamma+float(sv*np.transpose(sv))*1j)/vol/med
    if np.any(sv):
        return res-math.pi/gamma/gamma/vol-1/np.linalg.norm(sv)
    else:
        return res-math.pi/gamma/gamma/vol-2*gamma/math.sqrt(math.pi*la.det(dt))


def multi(a,b):
    if not a.shape[1] == b.shape[1]:
        return False
    else:
        c=[]
        for i in range(a.shape[1]):
            c.append(np.array(a[i]*b[0,i])[0])
        return np.mat(c)

def sampling(cell,prec):
    samples=[]
    for i in range(prec):
        i /= prec
        for j in range(prec):
            j /= prec
            for k in range(prec):
                k /= prec
                samples.append(cell[0]*i+cell[1]*j+cell[2]*k)
    return samples

def sampling_2d(cell,prec):
    samples=[]
    for i in range(prec+1):
        for j in range(prec-i+1):
            samples.append((cell[0]*i+cell[1]*j)/prec)
    return samples

def mono():
    size=2
    z_range=30
    for x in range(size)[1:]:
        for y in range(size)[x:]:
            for z in range(z_range)[28:]:
                vectors=multi(np.mat([[2.504651,0,0],[-1.252325,2.169091,0],[0,0,1]]),np.mat([x,y,z]))
                dielectric_tensor=np.mat([[24.945409/z+1,0,0],[0,24.945409/z+1,0],[0,0,z/(z-4.172816)]])
                vol=float(vectors[0]*np.transpose(np.cross(vectors[1],vectors[2])))
                rvectors=[]
                rvectors.append(2*math.pi*np.cross(vectors[1],vectors[2])/vol)
                rvectors.append(2*math.pi*np.cross(vectors[2],vectors[0])/vol)
                rvectors.append(2*math.pi*np.cross(vectors[0],vectors[1])/vol)
                current_gamma=1e10
                for gamma in range(40)[5:]:
                    gamma /= 20
                    current_radius=ewald_sum(gamma,1,vectors,vol,rvectors,dielectric_tensor,np.mat([[0,0,0]]))
                    flag=False
                    for radius in range(20)[1:]:
                        succeed_radius=ewald_sum(gamma,radius+1,vectors,vol,rvectors,dielectric_tensor,np.mat([[0,0,0]]))
                        if abs(current_radius-succeed_radius) < 1e-10:
                            succeed_gamma=current_radius
                            flag=True
                            print(x,y,z,gamma,radius,current_radius)
                            break
                        else:
                            current_radius=succeed_radius
                    if flag:
                        if abs(current_gamma-succeed_gamma) < 1e-10:
                            print(x,y,z,current_gamma)
                            break
                        else:
                            current_gamma=succeed_gamma
    return None


def iso():
    size=100
    dielectric_tensor=np.mat([[1,0,0],[0,1,0],[0,0,1]])
    for x in range(size)[1:2]:
        for y in range(size)[10:]:
            for z in range(size)[10:y+1]:
                y_resized=y/size
                z=z/size
                vectors=multi(np.mat([[1,0,0],[0,1,0],[0,0,1]]),np.mat([x,y_resized,z]))
                vol=float(vectors[0]*np.transpose(np.cross(vectors[1],vectors[2])))
                rvectors=[]
                rvectors.append(2*math.pi*np.cross(vectors[1],vectors[2])/vol)
                rvectors.append(2*math.pi*np.cross(vectors[2],vectors[0])/vol)
                rvectors.append(2*math.pi*np.cross(vectors[0],vectors[1])/vol)
                gamma_list=[]
                current_nei=1e10
                gamma_ref=1.8/np.power(vol,1/3)
                gamma_size=10
                for gamma in range(gamma_size):
                    gamma=((gamma/gamma_size)*0.4+0.8)*gamma_ref
                    current_radius=ewald_sum(gamma,1,vectors,vol,rvectors,dielectric_tensor,np.mat([[0,0,0]]))
                    flag=False
                    for radius in range(15)[1:]:
                        succeed_radius=ewald_sum(gamma,radius+1,vectors,vol,rvectors,dielectric_tensor,np.mat([[0,0,0]]))
                        # print(x,y_resized,z,gamma,radius+1,succeed_radius.real)
                        if abs(current_radius.real-succeed_radius.real) < 1e-5:
                            # print(x,y_resized,z,gamma,radius+1,succeed_radius)
                            gamma_list.append([gamma,radius+1,succeed_radius.real])
                            flag=True
                            break
                        else:
                            current_radius=succeed_radius
                    if flag:
                        if radius > current_nei:
                            break
                        else:
                            current_nei=radius
                gamma_max=0
                gamma_min=1e10
                for i in gamma_list:
                    if i[1] == current_nei+1:
                        if i[0] > gamma_max:
                            gamma_max=i[0]
                            esum_max=i[2]
                        if i[0] < gamma_min:
                            gamma_min=i[0]
                            esum_min=i[2]
                print(x,y_resized,z,vol,gamma_max,esum_max,gamma_min,esum_min,current_nei +1)
    # file_object.close()
    return None

def hbn():
    size=20
    #dielectric_tensor=np.mat([[4.746712,0,0],[0,4.746712,0],[0,0,2.679113]])
    dielectric_tensor=np.mat([[5.675642,0,0],[0,5.675642,0],[0,0,5.849838]])
    d_avg=math.pow(dielectric_tensor[0,0]*dielectric_tensor[1,1]*dielectric_tensor[2,2],1/3)
    for x in range(size)[1:]:
        for y in range(size)[x:]:
            #y=7+7*(y/size-1/2)
            for z in range(size)[1:]:
                #z=2+1.8*(z/size-1/2)
                #vectors=multi(np.mat([[2.504651,0,0],[-1.252325,2.169091,0],[0,0,6.657263]]),np.mat([x,y,z]))
                vectors=multi(np.mat([[3.1810592099274531,0,0],[-1.5905296050126307,2.7548780867520462,0],[0,0,5.1814313600938711]]),np.mat([x,y,z]))
                vol=float(vectors[0]*np.transpose(np.cross(vectors[1],vectors[2])))
                rvectors=[]
                rvectors.append(2*math.pi*np.cross(vectors[1],vectors[2])/vol)
                rvectors.append(2*math.pi*np.cross(vectors[2],vectors[0])/vol)
                rvectors.append(2*math.pi*np.cross(vectors[0],vectors[1])/vol)
                gamma_list=[]
                current_nei=1e10
                for gamma in range(200)[1:]:
                    gamma /= 20
                    current_radius=ewald_sum(gamma,1,vectors,vol,rvectors,dielectric_tensor,np.mat([[0,0,0]]))
                    flag=False
                    for radius in range(15)[1:]:
                        succeed_radius=ewald_sum(gamma,radius+1,vectors,vol,rvectors,dielectric_tensor,np.mat([[0,0,0]]))
                        if abs(current_radius-succeed_radius) < 1e-5:
                            #print(x,y,z,gamma,radius+1,succeed_radius)
                            gamma_list.append([gamma,radius+1,succeed_radius])
                            flag=True
                            break
                        else:
                            current_radius=succeed_radius
                    if flag:
                        if radius > current_nei:
                            break
                        else:
                            current_nei=radius
                gamma_max=0
                gamma_min=1e10
                for i in gamma_list:
                    if i[1] == current_nei+1:
                        if i[0] > gamma_max:
                            gamma_max=i[0]
                            esum_max=i[2]
                        if i[0] < gamma_min:
                            gamma_min=i[0]
                            esum_min=i[2]
                print(x,y,z,vol,d_avg,gamma_max,esum_max,gamma_min,esum_min,current_nei +1)
    return None

def pot():
    size=2
    sample_size=20
    dielectric_tensor=np.mat([[1,0,0],[0,1,0],[0,0,1]])
    d_avg=math.pow(dielectric_tensor[0,0]*dielectric_tensor[1,1]*dielectric_tensor[2,2],1/3)
    for x in range(size)[1:]:
        for y in range(size)[x:]:
            for z in range(size)[y:]:
                vectors=multi(np.mat([[1,0,0],[0,3,0],[0,0,5]]),np.mat([x,y,z]))
                vol=float(vectors[0]*np.transpose(np.cross(vectors[1],vectors[2])))
                rvectors=[]
                rvectors.append(2*math.pi*np.cross(vectors[1],vectors[2])/vol)
                rvectors.append(2*math.pi*np.cross(vectors[2],vectors[0])/vol)
                rvectors.append(2*math.pi*np.cross(vectors[0],vectors[1])/vol)
                sample_list=sampling_2d(np.mat([[0,0,0.5],[0,0.5,0.5]]),sample_size)
                for i in sample_list:
                    current_gamma=1e10
                    for gamma in range(20)[1:]:
                        gamma /= 10
                        current_radius=ewald_sum(gamma,1,vectors,vol,rvectors,dielectric_tensor,i)
                        flag=False
                        for radius in range(15)[1:]:
                            succeed_radius=ewald_sum(gamma,radius+1,vectors,vol,rvectors,dielectric_tensor,i)
                            # print(x,y,z,i[0,0],i[0,1],i[0,2],gamma,radius+1,succeed_radius.real,succeed_radius.imag)
                            if abs(current_radius.real-succeed_radius.real) < 1e-8:
                                succeed_gamma=current_radius
                                print(x,y,z,i[0,0],i[0,1],i[0,2],gamma,radius+1,succeed_radius.real,succeed_radius.imag)
                                flag=True
                                break
                            else:
                                current_radius=succeed_radius
                        if flag:
                            if abs(current_gamma.real-succeed_gamma.real) < 1e-5:
                                print(x,y,z,i[0,0],i[0,1],i[0,2],current_gamma.real,current_gamma.imag)
                                break
                            else:
                                current_gamma=succeed_gamma
    return None

print("14.414.414.414.414.414.414.4")
#hbn()
pot()


'''

for sc and iso
vectors=np.mat([[1,0,0],[0,1,0],[0,0,1]])
dielectric_tensor=np.mat([[1,0,0],[0,1,0],[0,0,1]])

for hBN
vectors=np.mat([[2.504651,0,0],[-1.252325,2.169091,0],[0,0,6.657263]])
dielectric_tensor=np.mat([[4.746712,0,0],[0,4.746712,0],[0,0,2.679113]])

for hBN monolayer
vectors=np.mat([[2.514292,0,0],[-1.257146,2.177441,0],[0,0,9.988486]])
dielectric_tensor=np.mat([[2.266586,0,0],[0,2.266586,0],[0,0,1.261810]])
ion_contribution=np.mat([[1.768008,0.000003,0],[0.000003,1.770692,0],[0,0,0.442994]])
'''

