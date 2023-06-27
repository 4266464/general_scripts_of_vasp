# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:18:39 2021

@author: dugue
"""

def wf_sum(path,band):
    import numpy as np
    text1=path+"\WF_SOC_SQUARED_B"+(4-len(str(band)))*"0"+str(band)+"_K0001_DW.vasp"
    text2=path+"\WF_SOC_SQUARED_B"+(4-len(str(band)))*"0"+str(band)+"_K0001_UP.vasp"
    sum_text=path+"\WF_SOC_SQUARED_B0"+str(band)+"_K0001.vasp"
    with open(text1,"r") as f:
        lines1=f.readlines()
    with open(text2,"r") as f:
        lines2=f.readlines()
    wf_start=False
    with open(sum_text,"w") as f:
        for i,j in zip(lines1,lines2):
            if len(i)<3 and wf_start==False:
                f.write(i)
                wf_start="ready"
            elif wf_start=="ready":
                f.write(i)
                wf_start=True
            elif wf_start==False:
                f.write(i)
            elif len(i)<3 and wf_start==True:
                f.write(i)
            else:
                for k in (np.array([float(k) for k in i.replace("\n","").split(" ") if k])+
                          np.array([float(k) for k in j.replace("\n","").split(" ") if k])):
                    f.write("   "+str(k))
                f.write("\n")

wf_sum(r"e:\Desktop",1305)