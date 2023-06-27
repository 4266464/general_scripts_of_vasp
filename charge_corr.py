#!/share/home/ckduan/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 19:55:38 2021

@author: dugue
"""
# the final correction should multiply 14.4
# run in queue
import numpy as np
from numpy import linalg as la
import scipy.special
import math
from old_poscar import poscar
import read_diele
import argparse
import pprint

def nei(radius, vector):
    list_nei = []
    for i in range(-radius, radius + 1):
        for j in range(-radius, radius + 1):
            for k in range(-radius, radius + 1):
                current_vector = np.array([i, j, k]).dot(vector)
                dis = la.norm(current_vector)
                if not dis > radius * min([la.norm(l) for l in vector]) and dis:
                    list_nei.append(current_vector)
    return list_nei


def ewald_sum(gamma, radius, vectors, dielectric_tensor, sv):
    res = 0
    vol = vectors[0].dot(np.cross(vectors[1], vectors[2]))
    rvectors = 2 * math.pi / vol * np.array([np.cross(vectors[1], vectors[2]),
                                             np.cross(vectors[2], vectors[0]),
                                             np.cross(vectors[0], vectors[1])])
    list_nei_real = nei(radius, vectors)
    if np.any(sv):
        list_nei_real.append(np.array([0, 0, 0]))
        for i in range(len(list_nei_real)):
            list_nei_real[i] = list_nei_real[i] - sv
    list_nei_rec = nei(radius, rvectors)
    for i in list_nei_real:
        med = math.sqrt(i.dot(la.inv(dielectric_tensor).dot(i)))
        res += scipy.special.erfc(gamma * med) / med / math.sqrt(la.det(dielectric_tensor))
    for i in list_nei_rec:
        med = i.dot(dielectric_tensor.dot(np.transpose(i)))
        res += 4 * math.pi * np.exp(-med / 4 / gamma / gamma + sv.dot(np.transpose(sv)) * 1j) / vol / med
    if np.any(sv):
        return res - math.pi / gamma / gamma / vol - 1 / np.linalg.norm(sv)
    else:
        return res - math.pi / gamma / gamma / vol - 2 * gamma / math.sqrt(math.pi * la.det(dielectric_tensor))


def iso(vectors, dielectric_tensor):
    vol = vectors[0].dot(np.cross(vectors[1], vectors[2]))
    gamma_list = []
    current_nei = 1e5
    gamma_ref = 1.8 / np.power(vol, 1 / 3)
    for i in np.linspace(0.8, 1.2, 11):
        gamma = i * gamma_ref
        current_radius = ewald_sum(gamma, 1, vectors, dielectric_tensor, np.array([0, 0, 0]))
        flag = False
        for radius in range(15)[1:]:
            succeed_radius = ewald_sum(gamma, radius + 1, vectors, dielectric_tensor, np.array([0, 0, 0]))
            if abs(current_radius.real - succeed_radius.real) < 1e-4:
                gamma_list.append([gamma, radius + 1, succeed_radius.real])
                flag = True
                break
            else:
                current_radius = succeed_radius
        if flag:
            if radius > current_nei:
                break
            else:
                current_nei = radius
    sum_list = list(zip(*gamma_list))[2]
    if max(sum_list) - min(sum_list) < 1e-4:
        return max(sum_list)
    else:
        return gamma_list


# def read_diele(text_path):
#     with open(text_path, "r") as f:
#         lines = f.readlines()
#     diele_ele = {}
#     diele_ion = {}
#     for i in range(len(lines)):
#         if "MACROSCOPIC STATIC DIELECTRIC TENSOR" in lines[i]:
#             if "IONIC" in lines[i]:
#                 diele_ion[lines[i].split("/")[0]] = np.array(
#                     [[float(k) for k in [k for k in lines[i + j + 2].split(" ") if k][1:4]] for j in range(3)])
#
#             elif not lines[i].split("/")[0] in diele_ele:
#                 diele_ele[lines[i].split("/")[0]] = np.array(
#                     [[float(k) for k in [k for k in lines[i + j + 2].split(" ") if k][1:4]] for j in range(3)])
#
#     print(diele_ele)
#     print(diele_ion)
#     return diele_ele, diele_ion


# def read_vector(text_path):
#     with open(text_path, "r") as f:
#         lines = f.readlines()
#     vector_list = []
#     tensor = False
#     for i in lines:
#         if "==>" in i:
#             vector_list.append([i.replace("==> ", "").split("/")[0]])
#         elif "1.0" in i:
#             tensor = []
#         elif type(tensor) == list:
#             if len(tensor) < 3:
#                 tensor.append([float(j) for j in i.split(" ") if j])
#
#                 if len(tensor) == 3:
#                     vector_list[-1].append(np.array(tensor))
#                     tensor = False
#     print(vector_list)
#     return vector_list

def iso_from_file(poscar_path, outcar_path):
    aposcar = poscar()
    aposcar.read(poscar_path)
    dielectric_tensor_dict = read_diele.read_dielectric_tensor_from_outcar(outcar_path)
    electric_contribution = dielectric_tensor_dict[
        'MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)']
    ionic_contribution = dielectric_tensor_dict[
        'MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION']
    pprint.pprint(f"cell paramaters:\n{aposcar.axis}")
    pprint.pprint(f"electric contribution:\n{electric_contribution}")
    pprint.pprint(f"ionic contribution:\n{ionic_contribution}")
    pprint.pprint(f"both:\n{electric_contribution + ionic_contribution}")
    return iso(aposcar.axis, electric_contribution + ionic_contribution) * 14.4 / 3


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--poscar_path", default="poscar", help="the path of **supercell** poscar, default 'poscar' ")
    parser.add_argument("-o", "--outcar_path", default="OUTCAR", help="the path of OUTCAR, default 'OUTCAR' ")
    args = parser.parse_args()
    print(iso_from_file(args.poscar_path, args.outcar_path))
