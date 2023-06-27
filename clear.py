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

def clear(del_list,del_path):
    for i,j,k in os.walk(del_path):
        for l in del_list:
            if l in k:
                linux_command("rm "+i+"/"+l)
        if not "relax" in i and "WAVECAR" in k:
            linux_command("rm "+i+"/WAVECAR")
    return None

clear(['PCDAT','EIGENVAL','DOSCAR','REPORT',
       'XDATCAR','vasprun.xml','IBZKPT','CHG','CHGCAR'],
      ".")