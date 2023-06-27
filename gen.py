# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 20:45:54 2021

@author: dugue
"""

def replace(text_path,replace_list):
    with open(text_path,"r") as f:
        lines=f.readlines()
    for i in replace_list:
        with open(i["path"],"w") as f:
            for j in lines:
                for k in i:
                    if not k=="path":
                        j=j.replace(k,i[k])
                f.write(j)

def text2list(text_path):
    with open(text_path,"r") as f:
        lines=f.readlines()
    key=[i for i in lines[0].replace("\n","").split(",") if not i[0]=="#" and i]
    replace_list=[]
    for i in lines[1:]:
        if not "#" in i:
            temp_dict={}
            value=[j for j in i.replace("\n","").split(",") if j]
            for j,k in zip(key,value):
                temp_dict[j]=k
            replace_list.append(temp_dict)
    return replace_list

replace(r"e:\Desktop\temlplate\Bi and Sb\DHP\temp_DHP.csv",
        text2list(r"e:\Desktop\template.csv"))
                    