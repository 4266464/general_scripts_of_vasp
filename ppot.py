# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 22:30:12 2021

@author: dugue
"""

def atom_pot(vhar,path):
    import re
    import matplotlib.pyplot as plt
    with open(f"{vhar}/PLANAR_AVERAGE.dat","r") as f:
        lines=f.readlines()
    match_result=re.search("#(.+), (.+)",lines[0])
    if not match_result:
        match_result=re.search("#(.+) (.+)",lines[0])
    xlabel,ylabel=match_result.groups()
    x=[xlabel]
    y=[ylabel]
    for j in lines[1:]:
        temp=[i for i in j.replace("\n","").split(" ") if i]
        x.append(float(temp[0]))
        y.append(-float(temp[1]))
    
    avg=sum(y[1:])/len(y[1:])
    print(y[1],avg)
    plt.plot(x[1:],y[1:])
    plt.plot([x[1],x[-1]],[y[1],y[-1]],"b--")
    plt.plot([x[1],x[-1]],[avg,avg],"r--")
    plt.plot([x[1]+0.1,x[1]+0.1],[y[1],avg],"k-^")
    plt.annotate(str(-round(y[1]-avg,3))+"eV",xy=(x[1],(y[1]-avg)/2),xytext=(10,0),textcoords='offset points',fontsize=12)
    plt.xlabel(x[0])
    plt.ylabel(y[0])
    #plt.show()
    plt.savefig(path+".png",dpi=300)
    plt.close('all')
    
    
    # count={}
    # for j in [round(j,3) for j in y[1:]]:
    #     if j in count:
    #         count[j]+=1
    #     else:
    #         count[j]=1
    
    # count=1
    # ysum=y[1]
    # temp_count=[]
    # temp_ysum=[]
    # for j in range(len(x))[2:-1]:
    #     count+=1
    #     ysum+=y[j]
    #     if (y[j]-y[j-1])/(x[j]-x[j-1])*(y[j+1]-y[j])/(x[j+1]-x[j])<0:
    #         temp_count.append(count)
    #         temp_ysum.append(ysum)
            
    #         count=0
    #         ysum=0
    
    
    # for j,k in zip(temp_count,temp_ysum):
    #     print(str(j)+" "+str(k))
    
    # good_range=3
    # good_list=[]
    # for j in range(len(x))[good_range+1:len(x)-good_range]:
    #     if abs((y[j]-y[j-good_range])/(x[j]-x[j-good_range])-
    #            (y[j]-y[j+good_range])/(x[j]-x[j+good_range]))<0.01:
    #         good_list.append([x[j],y[j]])
    # if len(good_list)>5:
    #     print("avg middle")
    #     print(str(sum(list(zip(*good_list))[1])/len(good_list))+" "+str(good_list[int(len(good_list)/2)][1]))
    return None

def slab_pot(vhar,path):
    import re
    import matplotlib.pyplot as plt
    with open(vhar,"r") as f:
        lines=f.readlines()
    match_result=re.search("#(.+), (.+)",lines[0])
    if not match_result:
        match_result=re.search("#(.+) (.+)",lines[0])
    xlabel,ylabel=match_result.groups()
    x=[xlabel]
    y=[ylabel]
    for j in lines[1:]:
        temp=[i for i in j.replace("\n","").split(" ") if i]
        x.append(float(temp[0]))
        y.append(float(temp[1]))
    
    #avg=sum(y[1:])/len(y[1:])
    plt.plot(x[1:],y[1:])
    # print(y[1],avg)
    # plt.plot([x[1],x[-1]],[y[1],y[-1]],"b--")
    # plt.plot([x[1],x[-1]],[avg,avg],"r--")
    # plt.plot([x[1]+0.1,x[1]+0.1],[y[1],avg],"k-^")
    # plt.annotate(str(round(y[1]-avg,3))+"eV",xy=(x[1],(y[1]-avg)/2),xytext=(10,0),textcoords='offset points',fontsize=12)
    plt.xlabel(x[0])
    plt.ylabel(y[0])
    #plt.show()
    #
    #plt.close('all')
    
    count=1
    ysum=y[1]
    temp_count=[]
    temp_ysum=[]
    slope=[(y[2]-y[1])/(x[2]-x[1])]
    for j in range(len(x))[2:-1]:
        count+=1
        ysum+=y[j]
        slope.append((y[j+1]-y[j])/(x[j+1]-x[j]))
        if slope[-1]*slope[-2]<-1e-8:
            if 25>count>10:
                temp_count.append(count)
                temp_ysum.append(ysum)
            count=0
            ysum=0
    #print(slope)
    count+=1
    ysum+=y[-1]
    if 25>count>10:
        temp_count.append(count)
        temp_ysum.append(ysum)
    print(temp_count,temp_ysum)
    
    avg=sum(temp_ysum)/sum(temp_count)
    plt.plot([x[1],x[-1]],[avg,avg],"r--")
    
    good_range=3
    good_list=[]
    for j in range(len(x))[good_range+1:len(x)-good_range]:
        if abs((y[j]-y[j-good_range])/(x[j]-x[j-good_range])-
                (y[j]-y[j+good_range])/(x[j]-x[j+good_range]))<0.01:
            good_list.append([x[j],y[j]])
    if len(good_list)>5:
        vac=good_list[int(len(good_list)/2)][1]
        plt.plot([x[1],x[-1]],[vac,vac],"r--")
        plt.plot([(x[1]+x[-1])/2,(x[1]+x[-1])/2],[vac,avg],"k-v")
        plt.annotate(str(-round(vac-avg,2))+"eV",xy=((x[1]+x[-1])/2,(vac+avg)/2),xytext=(10,0),textcoords='offset points',fontsize=12)
    
    plt.savefig(path+".png",dpi=300)
    print(vac,avg)
    #plt.show()
    plt.close('all')
    
    return None

if __name__ == '__main__':
    #atom_pot("cs.dat")
    slab_pot(r"e:\Desktop\10.12\PLANAR_AVERAGE (8).dat","hf3")