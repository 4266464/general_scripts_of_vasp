# 文件格式
# ion 1
#   spin1
#   spin2
# ion2
# ...

o_list=[3,0,3,0,3,0,3,0,
        3,0,3,0,3,0,3,0,
        3,0,3,0,3,0,3,0,
        3,0,3,0,3,0,3,0,
        3,0,3,0]

no_list=[1,
            1,1,1,1,
            1,1,1,1,
        1,
            1,1,1,1,
            1,1,1,1,
        1,
            1,1,1,1,
            1,1,1,1,
        1,
            1,1,1,1,
            1,1,1,1,]

sort=[27,9,28,10,29,11,30,12,
      31,0,18,1,19,2,20,3,
      21,4,22,5,23,6,24,7,
      25,8,26,13,32,14,33,15,
      34,16,35,17]




import re
import numpy as np
pattern=re.compile("([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)")
with open("vasprun.xml","r") as f:
    lines=f.readlines()

up=[]
down=[]
l=[]
x=[]
flag=False
aflag=False
for i in lines:
    if type(flag)==str:
        if "<r>" in i:
            l.append(sum([float(j) for j in [j for j in i.split(" ") if j][2:-1]]))
            if not aflag:
                x.append(float([j for j in i.split(" ") if j][1]))
        else:
            if "1" in flag:
                up.append(l)
                flag=True
            else:
                down.append(l)
                flag=False
            l=[]
            aflag=True
            
    else:
        if "<set comment=\"ion" in i:
            flag=True
        elif "<set comment=\"spin" in i and flag==True:
            flag=i



o_up=up[:54]

o_down=down[:54]
no_up=up[54:]

no_down=down[54:]

# o_up=list(zip(*o_up))[219]
# print(o_up)
# no_up=list(zip(*no_up))[219]
# print(no_up)
# o_down=list(zip(*o_down))[219]
# print(o_down)
# no_up=list(zip(*no_up))[219]
# print(no_up)
sum_up=[]
sum_down=[]

ci=0
cj=0
sum_up=[]
sum_down=[]

for i,j in zip(o_list,no_list):
    tsum_up=np.zeros(1500)
    tsum_down=np.zeros(1500)
    if i:
        for k in range(ci+i)[ci:]:
            tsum_up+=np.array(o_up[k])
            tsum_down+=np.array(o_down[k])
        ci+=i
    if j:
        for k in range(cj+j)[cj:]:
            tsum_up+=np.array(no_up[sort[k]])
            tsum_down+=np.array(no_down[sort[k]])
        cj+=j
    sum_up.append(tsum_up)
    sum_down.append(tsum_down)

# sum_up=[sum_up[sort[i]] for i in range(len(sum_up))]
# sum_down=[sum_down[sort[i]] for i in range(len(sum_down))]

with open("output","w") as f:
    for i in range(len(sum_up[0])):
        f.write(str(x[i])+" ")
        for j in range(len(sum_up)):
            f.write(str(sum_up[j][i])+" "+str(sum_down[j][i])+" ")
        f.write("\n")
