# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 02:14:04 2021

@author: dugue
"""

import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
import numpy as np
import json
import random
import os


def matproj_formation_energy(formula):
    prop = ['unit_cell_formula', 'material_id', 'energy_per_atom', 'formation_energy_per_atom']
    get_data = False
    with open("mpb", "a+") as f:
        f.seek(0)
        for i in f.readlines():
            if type(get_data) == list:
                if "search " in i:
                    break
                else:
                    get_data.append(json.loads(i.replace("\n", "")))
            elif i.replace("search ", "").replace("\n", "") == formula:
                get_data = []

    if type(get_data) == bool:
        get_data = []
        with open("mpb", "a") as f:
            with MPRester("VUswQqBEWe4VFBZD25") as m:
                f.write("search " + formula + "\n")
                for i in m.get_data(formula):
                    get_data.append(i)
                    f.write(json.dumps(i) + "\n")

    prop_dict = {}
    for i in get_data:
        agcd = gcd(list(i['unit_cell_formula'].values()))
        if agcd:
            for j in i['unit_cell_formula']:
                i['unit_cell_formula'][j] /= agcd
        item = tuple(sorted(list(i['unit_cell_formula'].items()), key=lambda x: x[0]))
        if not item in prop_dict or i['formation_energy_per_atom'] < prop_dict[item][-1]:
            prop_dict[item] = [i[j] for j in prop[1:]]

    return prop_dict


def gcd(alist):
    blist = list(set(alist[:]))
    clist = blist[:]
    while True:
        if 1 in blist:
            return False
        else:
            for i in blist:
                for j in blist:
                    sub = abs(i - j)
                    if sub and not sub in clist:
                        clist.append(sub)
            if len(clist) == len(blist):
                return min(blist)
            else:
                blist = clist[:]


def formula2tuple(formula):
    chem_dict = {}
    ele = ""
    num = 0.
    for i in formula:
        if i.isupper():
            if ele:
                if not num:
                    chem_dict[ele] = 1.
                else:
                    chem_dict[ele] = num
            ele = i
            num = 0.
        elif i.islower():
            ele += i
        elif i.isdigit:
            num = num * 10 + float(i)
        else:
            print("illegal " + i)

    if not num:
        chem_dict[ele] = 1.
    else:
        chem_dict[ele] = num
    chem_dict = tuple(sorted(list(chem_dict.items()), key=lambda x: x[0]))

    return chem_dict


def tuple2formula(chem_tuple):
    formula = ""
    for i in chem_tuple:
        for j in i:
            if type(j) == float:
                if j == 1:
                    pass
                else:
                    formula += str(int(j))
            else:
                formula += j
    return formula


def line2point(point1, point2):
    if abs(point1[0] - point2[0]) < 1e-5:
        temp = [(point1[1] - point2[1]) / 1e-5, (point1[0] * point2[1] - point2[0] * point1[1]) / 1e-5]
    else:
        temp = [(point1[1] - point2[1]) / (point1[0] - point2[0]),
                (point1[0] * point2[1] - point2[0] * point1[1]) / (point1[0] - point2[0])]
    return temp


abc = [[1 / 3 ** 0.5, 1], [0, 0], [2 / 3 ** 0.5, 0]]
labc = [line2point(*i) for i in [[j for j in abc if not j == i] for i in abc]]


def parallel_line(ref_line, distance):
    if ref_line == labc[0]:
        return [ref_line[0], ref_line[1] + distance * (ref_line[0] ** 2 + 1) ** 0.5]
    else:
        return [ref_line[0], ref_line[1] - distance * (ref_line[0] ** 2 + 1) ** 0.5]


def distancetoline(line, xy):
    temp = (xy[0] * line[0] + line[1] - xy[1]) / (line[0] ** 2 + 1) ** 0.5
    return abs(temp)


def crosspoint(line1, line2):
    if abs(line1[0] - line2[0]) < 1e-5:
        temp = [-(line1[1] - line2[1]) / 1e-5, (line1[0] * line2[1] - line2[0] * line1[1]) / 1e-5]
    else:
        temp = [-(line1[1] - line2[1]) / (line1[0] - line2[0]),
                (line1[0] * line2[1] - line2[0] * line1[1]) / (line1[0] - line2[0])]
    return temp


def tritoxy(tri):
    if type(tri[0]) == bool:
        xy = crosspoint(parallel_line(labc[1], tri[1]), parallel_line(labc[2], tri[2]))
    elif type(tri[1]) == bool:
        xy = crosspoint(parallel_line(labc[0], tri[0]), parallel_line(labc[2], tri[2]))
    else:
        xy = crosspoint(parallel_line(labc[0], tri[0]), parallel_line(labc[1], tri[1]))
    return xy


def xytotri(xy):
    temp = np.array([distancetoline(i, xy) for i in labc])
    temp /= np.linalg.norm(temp, ord=2)
    return temp


def readtri(tri, formula_tuple, formation_energy):
    return [i * j[1] / formation_energy if not type(i) == bool else i for i, j in zip(tri, formula_tuple)]


def printtri(tri, formula_tuple, formation_energy):
    return [i / j[1] * formation_energy for i, j in zip(tri, formula_tuple)]


def condition(chempot_list, formation_list):
    return sum([i * j for i, j in zip(chempot_list, formation_list)]) < formation_list[-1]


# SOC crystal and mocular
soc_corr = {"Bi": [-0.473, -0.808],
            "Cs": [-0.131, -0.143],
            "Na": [-0.000, -0.000],
            "Cl": [-0.007, -0.002]}

# energy_corr={"Cl":0.614,"Bi":-0.808,"O":0.687}
energy_corr = {}
#{"Cl": 0.614, "O": 0.687}

def chempot_sample(formula, fixed={}, start_point=[], other_element=[], soc=False, sample_times=100):
    file_name = formula
    # if fixed:
    #     file_name+="(fixed="+",".join([i+":"+str(fixed[i]) for i in fixed])+")"
    # if other_element:
    #     file_name+="(other_element="+",".join(other_element)+")"
    file_name += ".csv"
    with open(file_name, "w") as f:
        # "Cs2SnCl6" (('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0))
        formula_tuple = formula2tuple(formula)
        formation_mat = {}
        other_formation_mat = {}
        energy_per_atom = {}
        for i in range(2 ** (len(formula_tuple) + len(other_element)))[1:]:
            bicode = bin(i)
            bicode = bicode.replace("0b", "".join(
                ["0" for j in range(len(formula_tuple) + len(other_element) - len(bicode) + 2)]))
            ions = []
            for j, k in zip(bicode, list(dict(formula_tuple).keys()) + other_element):
                if j == "1":
                    ions.append(k)
            formation_dict = matproj_formation_energy("-".join(ions))
            # {'unit_cell_formula':['material_id','energy_per_atom','formation_energy_per_atom']
            # {(('Cl', 3.0), ('Cs', 1.0), ('Sn', 1.0)): ['mp-27394', -3.533219568, -1.8115027265862065], 
            # (('Cl', 5.0), ('Cs', 1.0), ('Sn', 2.0)): ['mp-30164', -3.551530498125, -1.6656268253663797],
            # (('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0)): ['mp-608555', -3.3898603444444446, -1.922349615651341]}

            for j in formation_dict:
                if len(j) == 1:
                    energy_per_atom[j[0][0]] = formation_dict[j][1] + (
                        energy_corr[j[0][0]] if j[0][0] in energy_corr else 0)
                else:
                    if j==(('Al', 4.0), ('Ca', 1.0), ('O', 7.0)):
                        continue
                    temp_formation_mat = []
                    formation_energy = sum(list(zip(*j))[1]) * formation_dict[j][-1]
                    if soc:
                        for k in dict(j):
                            if k in soc_corr:
                                formation_energy += dict(j)[k] * (soc_corr[k][0] - soc_corr[k][1])

                    for k in list(dict(formula_tuple).keys()) + other_element:
                        if k in dict(j):
                            if k in fixed:
                                formation_energy -= dict(j)[k] * fixed[k]
                            else:
                                temp_formation_mat.append(dict(j)[k])
                        else:
                            if k in fixed or k in other_element:
                                pass
                            else:
                                temp_formation_mat.append(0)

                    temp_formation_mat.append(formation_energy)
                    temp_formation_mat.append(tuple2formula(j))
                    inter = set(other_element).intersection(set(dict(j)))
                    if len(inter) == 1:
                        if not list(inter)[0] in other_formation_mat:
                            other_formation_mat[list(inter)[0]] = {}
                        other_formation_mat[list(inter)[0]][j] = temp_formation_mat
                    elif len(inter) == 0:
                        formation_mat[j] = temp_formation_mat

        print("element energy(soc&anion corr)")
        for i in energy_per_atom:
            print(f"{i:7} {energy_per_atom[i]:5.3f}")

        print("material", end="       ")
        real_formula_tuple = tuple([i for i in formula_tuple if not i[0] in fixed])
        for i in real_formula_tuple:
            print(i[0], end="  ")
        print("formation")

        for i in formation_mat:
            temp = f"{formation_mat[i][-1]:12}"
            for j in formation_mat[i][:-2]:
                temp += f"{j:4.0f}"
            temp += f"{formation_mat[i][-2]:8.2f}"
            print(temp)

        if other_formation_mat:
            print("material", end="       ")
            for i in real_formula_tuple:
                print(i[0], end="  ")
            for i in other_element:
                print(i, end="  ")
            print("formation")
            for i in other_formation_mat:
                for j in other_formation_mat[i]:
                    temp = f"{other_formation_mat[i][j][-1]:12}"
                    for k in other_formation_mat[i][j][:-2]:
                        temp += f"{k:4.0f}"
                    temp += f"{other_formation_mat[i][j][-2]:8.2f}"
                    print(temp)

        chempot_min = [formation_mat[formula_tuple][-2] / i for i in formation_mat[formula_tuple][:-2]]

        f.write(",".join([str(i[0]) for i in real_formula_tuple]) + "\n")
        f.write(",".join([str(i) for i in chempot_min]) + "\n")

        for i in formation_mat:
            f.write(",".join([str(j) for j in formation_mat[i]]) + "\n")

        if len(real_formula_tuple) == 3:
            line_mat = {}
            ax1 = plt.subplot(111, aspect=1)
            for i in formation_mat:
                if not i == formula_tuple and formation_mat[i][-2] < 0:
                    num_notzero = len([1 for j in formation_mat[i][:-2] if not j == 0])
                    if num_notzero == 1:
                        line_mat[i] = []
                        for j in range(3):
                            if not formation_mat[i][j] == 0:
                                line_mat[i].append(0)
                                ratio = formation_mat[formula_tuple][j] / formation_mat[i][j]
                            else:
                                line_mat[i].append(formation_mat[formula_tuple][j])
                        line_mat[i].append(formation_mat[formula_tuple][3] - ratio * formation_mat[i][3])
                        line_mat[i].append(formation_mat[i][4])
                    elif num_notzero == 2:
                        line_mat[i] = formation_mat[i][:]
                    elif num_notzero == 3:
                        for j in range(3):
                            ratio = formation_mat[formula_tuple][j] / formation_mat[i][j]
                            line_mat[i] = []
                            for k in range(3):
                                line_mat[i].append(formation_mat[formula_tuple][k] - ratio * formation_mat[i][k])
                            if len([1 for j in line_mat[i] if not j == 0]) == 2:
                                break
                        line_mat[i].append(formation_mat[formula_tuple][3] - ratio * formation_mat[i][3])
                        line_mat[i].append(formation_mat[i][4])

            lines = {}
            for i in line_mat:
                end_points = []
                for j in range(3):
                    if not line_mat[i][j] == 0:
                        end_points.append([])
                        for k in range(3):
                            if j == k:
                                end_points[-1].append(line_mat[i][3] / line_mat[i][k])
                            elif line_mat[i][k] == 0:
                                end_points[-1].append(False)
                            else:
                                end_points[-1].append(0)

                #print(end_points)
                for j in end_points:
                    temp_a=readtri(j, real_formula_tuple, formation_mat[formula_tuple][-2])
                    #print(temp_a)
                    for k in range(len(temp_a)):
                        if type(temp_a[k])==bool:
                            temp_a[k]=1-sum(temp_a)
                            break
                    print(*temp_a)



                end_points[0] = tritoxy(readtri(end_points[0], real_formula_tuple, formation_mat[formula_tuple][-2]))
                end_points[1] = tritoxy(readtri(end_points[1], real_formula_tuple, formation_mat[formula_tuple][-2]))
                lines[line_mat[i][4]] = line2point(end_points[0], end_points[1])
                print(line_mat[i][4])
                # for j in list(zip(*end_points)):
                #     print(*j)
                ax1.plot(*list(zip(*end_points)), label=line_mat[i][4])
                ax1.plot(*list(zip(*end_points)), label=line_mat[i][4])
                ax1.annotate(line_mat[i][4], xy=end_points[0],xytext=(0, 0), textcoords='offset points', fontsize=2)
                ax1.annotate(line_mat[i][4], xy=end_points[1], xytext=(0, 0), textcoords='offset points', fontsize=2)
            print(*readtri([-3.105, -1.405, 0], real_formula_tuple, formation_mat[formula_tuple][-2]))

            for i, j in zip(formula_tuple, labc):
                lines[i[0] + "_rich"] = j

            lines_set = []
            for i in lines:
                for j in lines:
                    if not i == j and not {i, j} in lines_set:
                        lines_set.append({i, j})

        if start_point:
            chempot_list = start_point
        else:
            for i in range(100000000):
                res_formation = formation_mat[formula_tuple][-2]
                chempot_list = []
                for j in real_formula_tuple:
                    ratio = random.random()
                    chempot_list.append(res_formation * ratio / j[1])
                    res_formation *= (1 - ratio)

                for j in formation_mat:
                    if not j == formula_tuple:
                        if not condition(chempot_list, formation_mat[j][:-1]):
                            chempot_list = []
                            break

                if chempot_list:
                    break

        if chempot_list:
            if len(real_formula_tuple) == 3:
                ax1.scatter(*tritoxy(readtri(chempot_list, real_formula_tuple, formation_mat[formula_tuple][-2])),
                            color="red")
            f.write("sample\n" + ",".join([str(i) for i in chempot_list]) + "\n")
            chempot_range = {}
            for i in real_formula_tuple:
                chempot_range[(i[0], "max")] = chempot_list
                chempot_range[(i[0], "min")] = chempot_list

            ratio_list = [i / j for i, j in zip(chempot_list, chempot_min)]
            count = 0
            while count < sample_times:
                temp_ratio_list = [i + 0.1 * (2 * random.random() - 1) for i in ratio_list[:-1]]
                temp_ratio_list.append(1 - sum(temp_ratio_list))
                for i in temp_ratio_list:
                    if i > 1 or i < 0:
                        temp_ratio_list = []
                        break
                    else:
                        chempot_list = [j * k for j, k in zip(temp_ratio_list, chempot_min)]

                if not temp_ratio_list:
                    continue

                for i in formation_mat:
                    if not i == formula_tuple:
                        if not condition(chempot_list, formation_mat[i][:-1]):
                            chempot_list = []
                            break

                if chempot_list:
                    count += 1
                    ratio_list = temp_ratio_list
                    f.write(",".join([str(i) for i in chempot_list]) + "\n")

                    for i in range(len(real_formula_tuple)):
                        chempot_range[(real_formula_tuple[i][0], "max")] = (chempot_list
                                                                            if chempot_list[i] > chempot_range[
                            (real_formula_tuple[i][0], "max")][i]
                                                                            else chempot_range[
                            (real_formula_tuple[i][0], "max")])
                        chempot_range[(real_formula_tuple[i][0], "min")] = (chempot_list
                                                                            if chempot_list[i] < chempot_range[
                            (real_formula_tuple[i][0], "min")][i]
                                                                            else chempot_range[
                            (real_formula_tuple[i][0], "min")])

            chempot_range[("Cl","mid")]=[(i+j)/2 for i,j in zip(chempot_range[("Cl","max")],chempot_range[("Cl","min")])]

            print(chempot_range)
            result_dict = {}
            for i in chempot_range:
                if tuple(chempot_range[i]) in result_dict:
                    result_dict[tuple(chempot_range[i])]["tag"].append("_".join(i))
                else:
                    result_dict[tuple(chempot_range[i])] = {}
                    result_dict[tuple(chempot_range[i])]["tag"] = ["_".join(i)]
                    result_dict[tuple(chempot_range[i])]["chempot"] = {j[0]: k + energy_per_atom[j[0]] for j, k in
                                                                       zip(real_formula_tuple, chempot_range[i])}
            if other_element:
                for i in result_dict:
                    for j in other_formation_mat:
                        other_min = []
                        for k in other_formation_mat[j].values():
                            other_min.append([k[-1], (k[-2] - sum([i[l] * k[l] for l in range(len(i))])) / k[-3]])
                        other_min.sort(key=lambda x: x[-1])
                        result_dict[i]["tag"].append(j + "_" + other_min[0][0])
                        result_dict[i]["chempot"][j] = other_min[0][1] + energy_per_atom[j]
            if fixed:
                for i in result_dict:
                    for j in fixed:
                        result_dict[i]["chempot"][j] = fixed[j] + energy_per_atom[j]
            for i in result_dict:
                print(" ".join(result_dict[i]["tag"]))
                for j in result_dict[i]["chempot"]:
                    print(f"{j:5}{result_dict[i]['chempot'][j]:.3f}")
            if len(real_formula_tuple) == 3:
                lines = []
                for j, k in zip(real_formula_tuple, abc):
                    ax1.annotate(j[0] + "(" + ",".join([str(l) for l in printtri(xytotri(k), real_formula_tuple,
                                                                                 formation_mat[formula_tuple][
                                                                                     -2])]) + ")",
                                 k, xytext=(0, 0), textcoords='offset points', fontsize=8)
                    lines.append(k)

                for j in range(3):
                    for k in range(3)[j + 1:]:
                        temp = list(zip(lines[j], lines[k]))
                        print(*temp[0])
                        print(*temp[1])
                        ax1.plot(temp[0], temp[1], color="k")

                ax1.legend(fontsize=6)
                #plt.show()
                plt.savefig(formula + "_".join([i + "_" + str(fixed[i]) for i in fixed]) + ".png", dpi=500)
                plt.close('all')

            with open("js_" + formula, "w") as g:
                for i in result_dict:
                    g.write(json.dumps(result_dict[i]) + "\n")

        else:
            print("sample failed")

    return result_dict.values()


#chempot_sample("Cs2ZrCl6",sample_times=1000)
chempot_sample("Cs3Bi2Cl9", sample_times=1000)
#chempot_sample("Cs2HfCl6", sample_times=10000, other_element=["Sb", "Bi", "Se", "Te", "Zr"])

# sample=False
# with open("Cs2NaBiCl6.csv","r") as f:
#     for i in f.readlines():
#         if "range" in i:
#             sample.sort(key=lambda x:x[1])
#             mid=len(sample)//2
#             for j in [0,mid,-1]:
#             #for j in [0]:
#                 chempot_sample("Cs2NaBiCl6",
#                                 fixed={"Cl":sample[j][1]},
#                                 start_point=[sample[j][0],sample[j][2],sample[j][3]],
#                                 other_element=["O"],
#                                 soc=True)
#             break


#         if type(sample)==list:
#             sample.append([float(j) for j in i.replace("\n","").split(",")])

#         if "sample" in i:
#             sample=[]
