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
import re
from pymatgen.ext.matproj import MPRester
import json
import random


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


def formula2list(formula):
    element_list = []
    number_list = []
    re_formula = ""
    for i in formula:
        if i.isupper():
            re_formula += r"([A-Z][a-z]?)(\d*)"
    for i in re.match(re_formula, formula).groups():
        if i == "":
            number_list.append(1)
        elif i[0].isdigit():
            number_list.append(int(i))
        else:
            element_list.append(i)
    # if order:
    #     sorted_number_list=[]
    #     for i in order:
    #         if i in element_list:
    #             sorted_number_list.append(number_list[element_list.index(i)])
    #         else:
    #             sorted_number_list.append(0)
    #     return sorted_number_list
    # else:
    return element_list, number_list


def list2formula(number_list, element_list):
    formula = ""
    for i, j in zip(number_list, element_list):
        if i == 0:
            pass
        elif i == 1:
            formula += j
        else:
            formula += j + str(int(i))
    return formula


def p2l(point1, point2):
    if abs(point1[0] - point2[0]) < 1e-5:
        temp = [(point1[1] - point2[1]) / 1e-5, (point1[0] * point2[1] - point2[0] * point1[1]) / 1e-5]
    else:
        temp = [(point1[1] - point2[1]) / (point1[0] - point2[0]),
                (point1[0] * point2[1] - point2[0] * point1[1]) / (point1[0] - point2[0])]
    return temp


def l2p(line1, line2):
    if abs(line1[0] - line2[0]) < 1e-5:
        temp = [-(line1[1] - line2[1]) / 1e-5, (line1[0] * line2[1] - line2[0] * line1[1]) / 1e-5]
    else:
        temp = [-(line1[1] - line2[1]) / (line1[0] - line2[0]),
                (line1[0] * line2[1] - line2[0] * line1[1]) / (line1[0] - line2[0])]
    return temp


# abc = [[1 / 3 ** 0.5, 1], [0, 0], [2 / 3 ** 0.5, 0]]
# labc = [line2point(*i) for i in [[j for j in abc if not j == i] for i in abc]]
# def parallel_line(ref_line, distance):
#     if ref_line == labc[0]:
#         return [ref_line[0], ref_line[1] + distance * (ref_line[0] ** 2 + 1) ** 0.5]
#     else:
#         return [ref_line[0], ref_line[1] - distance * (ref_line[0] ** 2 + 1) ** 0.5]
def distance(point, line):
    return (point[1] - point[0] * line[0] - line[1]) / (line[0] ** 2 + 1) ** 0.5


# def tritoxy(tri):
#     if type(tri[0]) == bool:
#         xy = crosspoint(parallel_line(labc[1], tri[1]), parallel_line(labc[2], tri[2]))
#     elif type(tri[1]) == bool:
#         xy = crosspoint(parallel_line(labc[0], tri[0]), parallel_line(labc[2], tri[2]))
#     else:
#         xy = crosspoint(parallel_line(labc[0], tri[0]), parallel_line(labc[1], tri[1]))
#     return xy
# def xytotri(xy):
#     temp = np.array([distancetoline(i, xy) for i in labc])
#     temp /= np.linalg.norm(temp, ord=2)
#     return temp
# def readtri(tri, formula_tuple, formation_energy):
#     return [i * j[1] / formation_energy if not type(i) == bool else i for i, j in zip(tri, formula_tuple)]
# def printtri(tri, formula_tuple, formation_energy):
#     return [i / j[1] * formation_energy for i, j in zip(tri, formula_tuple)]
# def condition(chempot_list, formation_list):
#     return sum([i * j for i, j in zip(chempot_list, formation_list)]) < formation_list[-1]
# SOC crystal and mocular
soc_corr = {"Bi": [-0.473, -0.808],
            "Cs": [-0.131, -0.143],
            "Cl": [-0.007, -0.002]}
corr = {"Cl": 0.614, "Bi": -0.808, "O": 0.687, "F": 0.462}


def chempot_plot(formula, fixed={}, exclude=[]):
    order, background = formula2list(formula)
    # absolute chemical potential
    atom_energy = {}
    # the big matrix
    mat = []
    for i in range(2 ** len(order))[1:]:
        bicode = bin(i)
        bicode = bicode.replace("0b", (len(order) - len(bicode) + 2) * "0")
        element_combination = []
        for j, k in zip(bicode, order):
            if j == "1":
                element_combination.append(k)
        formation_dict = matproj_formation_energy("-".join(element_combination))

        # {'unit_cell_formula':['material_id','energy_per_atom','formation_energy_per_atom']
        # {(('Cl', 3.0), ('Cs', 1.0), ('Sn', 1.0)): ['mp-27394', -3.533219568, -1.8115027265862065],
        # (('Cl', 5.0), ('Cs', 1.0), ('Sn', 2.0)): ['mp-30164', -3.551530498125, -1.6656268253663797]}
        for j in formation_dict:
            # 单质
            if len(j) == 1:
                atom_energy[j[0][0]] = formation_dict[j][1] + (corr[j[0][0]] if j[0][0] in corr else 0)
            # 化合物
            else:
                _mat = []
                formation_energy = sum(list(zip(*j))[1]) * formation_dict[j][-1]
                for k in order:
                    _mat.append(dict(j)[k] if k in dict(j) else 0)
                for k in dict(j):
                    if k in soc_corr:
                        formation_energy += dict(j)[k] * (soc_corr[k][0] - soc_corr[k][1])
                _mat.append(formation_energy)
                mat.append(_mat)

    # print main energy data

    print("material", end="       ")
    for i in order:
        print(i, end="  ")
    print("formation")

    tex = ""
    for i in mat:
        if i[:-1] == background:
            background.append(i[-1])
            # continue
        tex += f"{list2formula(i[:-1], order):12}"
        for j in i[:-1]:
            tex += f"{j:4.0f}"
        tex += f"{i[-1]:8.2f}\n"
    print(tex)

    print("atom energy")
    for i in atom_energy:
        print(f"{i:7} {atom_energy[i]:5.3f}")
    print()

    print(mat, background)

    mat.remove(background)
    for i in mat:
        i.append(list2formula(i[:-1], order))

    # 删除固定元素
    if fixed:
        for i in fixed:
            index = order.index(i)
            del order[index]
            background[-1] -= background[index] * fixed[i]
            del background[index]
            for j in mat:
                j[-2] -= j[index] * fixed[i]
                del j[index]

    points = {}
    line_background = [-background[1] / background[2], background[3] / background[2]]

    points[tuple(l2p([0, 0], [1e5, 0]))] = [f"{order[1]}=0", f"{order[2]}=0"]
    points[tuple(l2p([0, 0], line_background))] = [f"{order[0]}=0", f"{order[2]}=0"]
    points[tuple(l2p(line_background, [1e5, 0]))] = [f"{order[0]}=0", f"{order[1]}=0"]

    # 映射到正三角形
    _mat = []
    for i in points:
        if not i == (0, 0):
            _mat.append(i)
    # trans_mat=np.linalg.inv(np.array(_mat)).dot(np.array([[0.5,3**0.5*0.5],[1,0]]))
    trans_mat = np.array([[1, 0], [0, 1]])

    lines = {}
    lines[f"{order[0]}=0"] = line_background
    lines[f"{order[1]}=0"] = [1e5, 0]
    lines[f"{order[2]}=0"] = [0, 0]

    frame = {}
    frame[f"{order[0]}=0"] = [l2p([0, 0], line_background), l2p(line_background, [1e5, 0])]
    frame[f"{order[1]}=0"] = [l2p([0, 0], [1e5, 0]), l2p([1e5, 0], line_background)]
    frame[f"{order[2]}=0"] = [l2p([1e5, 0], [0, 0]), l2p([0, 0], line_background)]

    for i in mat:
        ################
        # 不考虑特定化合物
        ################
        if i[-1] in exclude:
            continue
        # 消元
        elimi = [j - k / background[0] * i[0] for j, k in zip(i[1:-1], background[1:])]
        # 不满足条件的端点
        bad_points = []
        for j in points:
            if j[0] * elimi[0] + j[1] * elimi[1] > elimi[2]:
                bad_points.append(j)

        # 忽略更稳定的化合物
        if len(bad_points) == len(points):
            print(f"more stable compound found: {i[-1]}")
            continue
            # return False

        # tobecut=[]
        # for j in bad_points:
        #     #判断那些点和j的连线在图形上
        #     for k in points:
        #         if not j==k:
        #             _line=p2l(j,k)
        #             _distance=[distance(l,_line) for l in points if not l in [j,k]]
        #             if max(_distance)*min(_distance)>0:
        #                 tobecut.append([_line,j,k])

        if bad_points:
            # 转化当前条件为kb
            if elimi[1] == 0:
                elimi[1] = 1e-5
            elimi = [-elimi[0] / elimi[1], elimi[2] / elimi[1]]
            for j in bad_points:
                for k in points[j]:
                    for l in points:
                        # j 是要切掉的点
                        # k 是可能要切断的线, 对于二维情形有两条线
                        # l 是通过k和j相连的点
                        if k in points[l] and l not in bad_points:
                            # 求交点
                            points[tuple(l2p(lines[k], elimi))] = [i[-1], k]
                            break
                del points[j]
            lines[i[-1]] = elimi

        # if tobecut:
        #     # 转化当前条件为kb
        #     elimi=[-elimi[0]/elimi[1],elimi[2]/elimi[1]]
        #     new_points= {}
        #     for j in tobecut:
        #         # 求交点
        #         _point=l2p(j[0],elimi)
        #         # 判断交点是否在图形上
        #         if 0<(_point[0]-j[1][0])/(j[2][0]-j[1][0])<1:
        #             new_points[tuple(_point)]=i[-1]
        #     # 替代点
        #     for j in bad_points:
        #         del points[j]
        #
        #     points.update(new_points)
        #
        #     # 可能要画的线
        #     for j in new_points:
        #         lines.append([i[-1],elimi])

    for i in points:
        plt.scatter(*np.array(i).dot(trans_mat), s=30)
        # plt.annotate(np.array(i).dot(trans_mat),xy=np.array(i).dot(trans_mat))

    # 单质线
    for i in frame:
        _line = list(zip(np.array(frame[i][0]).dot(trans_mat), np.array(frame[i][1]).dot(trans_mat)))
        plt.plot(_line[0], _line[1], label=i, linewidth=3)

    # 构成最终区域的线
    # tobeplotted=[]
    # for i in points:
    #     for j in lines[3:]:
    #         if points[i]==j[0]:
    #             if not j[0] in tobeplotted:
    #                 tobeplotted[j[0]]=j[1]

    tobeplotted = []
    for i in points:
        for j in points[i]:
            if j not in tobeplotted and "=0" not in j:
                tobeplotted.append(j)

    # tobeplotted.append("Bi2O3")

    # Origin output start
    tex = ""
    for i in order:
        tex += f"μ\-({i})\t"
    tex += "\n"
    for i in order:
        tex += "eV\t"
    tex += "\n"
    print(tex)

    chempot_min = [background[-1] / i for i in background[:-1]]

    # print(background,chempot_min)
    # Origin output end

    def lower_number(matched):
        return f"\-({matched.group('number')})"

    for i in tobeplotted:
        _points = []
        for j in frame:
            _point = l2p(p2l(frame[j][0], frame[j][1]), lines[i])
            # 判断交点是否在图形上
            if 0 < (_point[0] - frame[j][0][0]) / (frame[j][1][0] - frame[j][0][0]) < 1:
                _points.append(_point)
        # print(_points,i)
        tex = ""
        for j in _points:
            j.insert(0, (background[-1] - sum([k * l for k, l in zip(j, background[1:])])) / background[0])
            # print(j)
            for k, l in zip(j, chempot_min):
                tex += f"{k / l} "
            del j[0]
            tex += re.sub("(?P<number>\d+)", lower_number, i) + "\n"
        tex += "\n"
        print(tex)

        print(_points, trans_mat)

        _line = list(zip(np.array(_points[0]).dot(trans_mat), np.array(_points[1]).dot(trans_mat)))
        plt.plot(_line[0], _line[1], label=i, linewidth=3)

    tex = ""
    for i in points:
        j = list(i)
        j.insert(0, (background[-1] - sum([k * l for k, l in zip(j, background[1:])])) / background[0])
        print(" ".join([str(k) for k in j]))
        for k, l in zip(j, chempot_min):
            tex += f"{k / l} "
        tex += f"{points[i][0]} {points[i][1]}\n"
    print(tex)

    # print(points)
    # plt.xlim(-10,0)

    plt.legend(fontsize=20)
    # plt.lattice([0.5, -19, 1, -19.5])
    # plt.show()
    plt.savefig(f"{formula}.png")
    plt.close("all")
    return True


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
            elif set(i.replace("search ", "").replace("\n", "").split("-")) == set(formula.split("-")):
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


def can_form(point, formation_list):
    return sum([i * j for i, j in zip(point, formation_list)]) - formation_list[-1]


# SOC crystal and mocular
soc_corr = {"Bi": [-0.473, -0.808],
            "Cs": [-0.131, -0.143],
            "Cl": [-0.007, -0.002]}
corr = {"Cl": 0.614, "Bi": -0.808, "O": 0.687, "F": 0.462}


def chempot_sample(formula, start_point=[], fixed={}, other_element=[], exclude=[], soc=False, sample_times=1000):
    # "Cs2SnCl6" (('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0))

    formula_tuple = formula2tuple(formula)
    formation_mat = {}
    other_formation_mat = {}
    atom_energy = {}
    exclude = [set(formula2tuple(k)) for k in exclude]
    for i in range(2 ** (len(formula_tuple) + len(other_element)))[1:]:
        bicode = bin(i)
        bicode = bicode.replace("0b", (len(formula_tuple) + len(other_element) - len(bicode) + 2) * "0")
        element_combination = []
        for j, k in zip(bicode, list(dict(formula_tuple).keys()) + other_element):
            if j == "1":
                element_combination.append(k)
        formation_dict = matproj_formation_energy("-".join(element_combination))
        # {'unit_cell_formula':['material_id','atom_energy','formation_atom_energy']
        # {(('Cl', 3.0), ('Cs', 1.0), ('Sn', 1.0)): ['mp-27394', -3.533219568, -1.8115027265862065],
        # (('Cl', 5.0), ('Cs', 1.0), ('Sn', 2.0)): ['mp-30164', -3.551530498125, -1.6656268253663797],
        # (('Cl', 6.0), ('Cs', 2.0), ('Sn', 1.0)): ['mp-608555', -3.3898603444444446, -1.922349615651341]}

        for j in formation_dict:
            # 单质
            if len(j) == 1:
                atom_energy[j[0][0]] = formation_dict[j][1] + (corr[j[0][0]] if j[0][0] in corr else 0)
            # 化合物
            else:
                # 去掉更稳定化合物
                if set(j) in exclude:
                    pass
                else:
                    temp_formation_mat = []
                    formation_energy = sum(list(zip(*j))[1]) * formation_dict[j][-1]
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
                    # 判断化合物是否含有杂质
                    inter = set(other_element).intersection(set(dict(j)))
                    # 只考虑包含一种杂质的化合物
                    if len(inter) == 1:
                        if not list(inter)[0] in other_formation_mat:
                            other_formation_mat[list(inter)[0]] = {}
                        other_formation_mat[list(inter)[0]][j] = temp_formation_mat
                    # 无杂质
                    elif len(inter) == 0:
                        formation_mat[j] = temp_formation_mat

    print("material", end="       ")
    # 去掉fixed
    real_formula_tuple = tuple([i for i in formula_tuple if not i[0] in fixed])

    for i in real_formula_tuple:
        print(i[0], end="  ")
    print("formation")

    ############################################
    ############################################
    # formation
    ############################################
    ############################################

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

    print("atom energy")
    for i in atom_energy:
        print(f"{i:7} {atom_energy[i]:5.3f}")

    # 各组分化学势下限
    chempot_min = [formation_mat[formula_tuple][-2] / i for i in formation_mat[formula_tuple][:-2]]
    # f.write(",".join([str(i[0]) for i in real_formula_tuple]) + "\n")
    # f.write(",".join([str(i) for i in chempot_min]) + "\n")

    # for i in formation_mat:
    #    f.write(",".join([str(j) for j in formation_mat[i]]) + "\n")

    if start_point:
        pass
    else:
        # 寻找抽样初始点
        for i in range(10000000):
            res_formation = formation_mat[formula_tuple][-2]
            start_point = []
            for j in real_formula_tuple:
                ratio = random.random()
                start_point.append(res_formation * ratio / j[1])
                res_formation *= (1 - ratio)

            for j in formation_mat:
                if not j == formula_tuple:
                    if can_form(start_point, formation_mat[j][:-1]) > 0:
                        start_point = []
                        break

            if start_point:
                break

    if start_point:
        print(f"start_point:{start_point}")
        chempot_range = {}
        for i in real_formula_tuple:
            chempot_range[(i[0], "rich")] = start_point
            chempot_range[(i[0], "poor")] = start_point

        ratio_list = [i / j for i, j in zip(start_point, chempot_min)]
        count = 0
        while count < sample_times:
            _ratio_list = [i + 0.05 * random.random() * (2 * random.random() - 1) for i in ratio_list[:-1]]
            _ratio_list.append(1 - sum(_ratio_list))
            for i in _ratio_list:
                if i > 1 or i < 0:
                    _ratio_list = []
                    break
                else:
                    point = [j * k for j, k in zip(_ratio_list, chempot_min)]

            if _ratio_list == []:
                continue

            for i in formation_mat:
                if not i == formula_tuple:
                    if can_form(point, formation_mat[i][:-1]) > 0:
                        point = []
                        break

            if point:
                count += 1
                # 打印进度
                # if count > 100-1 and count % 100 == 0:
                #     print(f"{count / sample_times * 100:.2f}%")
                ratio_list = _ratio_list
                # f.write(",".join([str(i) for i in start_point]) + "\n")

                for i in range(len(real_formula_tuple)):
                    chempot_range[(real_formula_tuple[i][0], "rich")] = (
                        point if point[i] > chempot_range[(real_formula_tuple[i][0], "rich")][i] else
                        chempot_range[(real_formula_tuple[i][0], "rich")])
                    chempot_range[(real_formula_tuple[i][0], "poor")] = (
                        point if point[i] < chempot_range[(real_formula_tuple[i][0], "poor")][i] else
                        chempot_range[(real_formula_tuple[i][0], "poor")])

        for i in real_formula_tuple:
            chempot_range[(i[0], "moderate")] = [(j + k) / 2 for j, k in
                                                 zip(chempot_range[(i[0], "rich")], chempot_range[(i[0], "poor")])]
        _average = []
        for i in chempot_range:
            if i[1] in ["rich", "poor"]:
                _average.append(chempot_range[i])
        chempot_range[("all", "average")] = [sum(i) / len(_average) for i in list(zip(*_average))]

        print(chempot_range)

        # 合并tag
        result_dict = {}
        for i in chempot_range:
            if tuple(chempot_range[i]) in result_dict:
                result_dict[tuple(chempot_range[i])]["tag"].append("_".join(i))
            else:
                result_dict[tuple(chempot_range[i])] = {}
                result_dict[tuple(chempot_range[i])]["tag"] = ["_".join(i)]
                result_dict[tuple(chempot_range[i])]["chempot"] = {j[0]: k + atom_energy[j[0]] for j, k in
                                                                   zip(real_formula_tuple, chempot_range[i])}
        if other_element:
            for i in result_dict:
                for j in other_formation_mat:
                    other_min = []
                    for k in other_formation_mat[j].values():
                        other_min.append([k[-1], (k[-2] - sum([i[l] * k[l] for l in range(len(i))])) / k[-3]])
                    other_min.sort(key=lambda x: x[-1])
                    result_dict[i]["tag"].append(j + "_" + other_min[0][0])
                    result_dict[i]["chempot"][j] = other_min[0][1] + atom_energy[j]

        if fixed:
            for i in result_dict:
                for j in fixed:
                    result_dict[i]["chempot"][j] = fixed[j] + atom_energy[j]

        for i in result_dict:
            print(" ".join(result_dict[i]["tag"]))
            for j in result_dict[i]["chempot"]:
                print(f"{j:3}{result_dict[i]['chempot'][j]:8.3f}")

        with open("js_" + formula, "a+") as f:
            for i in result_dict:
                f.write(json.dumps(result_dict[i]) + "\n")

    else:
        print("sample failed")

    return result_dict.values()


if __name__ == "__main__":
    # for i in np.linspace(-1.976,-1.966,10):
    #     if chempot_plot("Cs2NaBiCl6",fixed={"Cl":i}):
    #         print(i)
    #         break

    # chempot_plot("Cs2NaBiCl6",fixed={"Cl":-0.430},other_element=["O"])
    # chempot_plot("Cs2NaBiCl6",fixed={"Cl":-1.969})
    # chempot_plot("Cs2NaBiCl6",fixed={"Cl":-1.200})
    # chempot_plot("La3Ga5SnO14", fixed={"O": -1.5})
    # chempot_sample("YPO4")
    # chempot_plot("CaAl12O19")
    # chempot_plot("SrAl12O19")
    # chempot_plot("BiOCl")
    # chempot_plot("Cs2HfCl6")
    # chempot_plot("O3AlY")
    # chempot_plot("Y3Al5O12")
    # chempot_plot("Mg2Al4Si5O18", fixed={"O":-4.9})
    # chempot_plot("RbMgF3")
    # chempot_plot("Cs2ZnCl4")
    chempot_plot("Ca5Ga6O14")

    # for i in ["SrGa4O7", "BaGa4O7", "CaGa4O7", "Sr2P2O7", "CaAl4O7", "BaAl4O7", "SrAl4O7"]:
    # # for i in ["SrGa4O7", "BaGa4O7", "CaGa4O7", "Sr2P2O7", "CaAl4O7", "SrAl4O7"]:
    #     chempot_plot(i)
    #     # chempot_sample(i, sample_times=10000, other_element=["Cr"])
    # chempot_sample("BaAl4O7",exclude=["Ba7Al64O103"], sample_times=10000, other_element=["Cr"])
    chempot_sample("Ca5Ga6O14", sample_times=100000, exclude=["Ca3Ga4O9"], other_element=["Bi"])
