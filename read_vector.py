# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 16:29:39 2021

@author: dugue
"""

def read_vector(text_path):
    import numpy as np
    with open(text_path,"r") as f:
        lines=f.readlines()
    vector_list=[]
    tensor=False
    for i in lines:
        if "==>" in i:
            vector_list.append([i.replace("==> ","").split("/")[0]])
        elif "1.00000000000000" in i:
            tensor=[]
        elif type(tensor)==list:
            if len(tensor)<3:
                tensor.append([float(j) for j in i.split(" ") if j])
                if len(tensor)==3:
                    vector_list[-1].append(np.array(tensor))
                    tensor=False
                
    for i in vector_list:
        print(i[0],i[1][0][0])
    return None

read_vector("vector")
