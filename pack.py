# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:30:23 2021

@author: dugue
"""

import os

def linux_command(command):
    import os
    try:
        with os.popen(command, "r") as p:
            command_return=p.readlines()
    except:
        print(command)
    return command_return

#del_list=['PCDAT','EIGENVAL','DOSCAR','REPORT','XDATCAR','vasprun.xml','IBZKPT','CHG','CHGCAR']
def package(pack_list,pack_path):
    for i,j,k in os.walk(pack_path):
        if not "pack_" in i:
            log_list=[]
            if set(k).intersection(set(pack_list)):
                for l in k:
                    if ".log" in l and l.replace(".log","").isdigit():
                        log_list.append(int(l.replace(".log","")))
            if log_list:
                linux_command("mkdir pack_"+i)
                log_list.sort()
                linux_command("cp "+i+"/"+str(log_list[-1])+".log"
                              +" pack_"+i+"/"+str(log_list[-1])+".log")
                for l in k:
                    if l in pack_list:
                        linux_command("cp "+i+"/"+l+" pack_"+i+"/"+l)
            elif j:
                linux_command("mkdir pack_"+i)
    while True:
        flag=True
        for i,j,k in os.walk(pack_path):
            if not j and not k:
                linux_command("rm -r "+i)
                flag=False
        if flag:
            break
    return None

package(['CONTCAR','INCAR','OUTCAR','KPOINTS','POSCAR','POTCAR','jobid','PROCAR'],".")
