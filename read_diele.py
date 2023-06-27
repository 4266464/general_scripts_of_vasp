# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 13:46:54 2021

@author: dugue
"""

import numpy as np
import argparse
import pprint

def read_dielectric_tensor_from_outcar(outcar_path):
    dielectric_tensor_dict = {}

    with open(outcar_path, "r") as outcar:
        while True:
            line = outcar.readline()
            if not line:
                break
            if "DIELECTRIC TENSOR" in line:
                dielectric_tensor = []
                outcar.readline()
                for i in range(3):
                    dielectric_tensor.append(
                        [float(dielectric_component) for dielectric_component in outcar.readline().strip().split(" ") if
                         dielectric_component])
                dielectric_tensor_dict[line.strip()]=np.array(dielectric_tensor)
    return dielectric_tensor_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("outcar_path",help="the path of OUTCAR")
    args = parser.parse_args()
    pprint.pprint(read_dielectric_tensor_from_outcar(args.outcar_path))