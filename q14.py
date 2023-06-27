# -*- coding: utf-8 -*-
"""
Created on Mon May 31 21:38:56 2021

@author: dugue

version=0.1
"""

# !/home/phys/qif/anaconda3/bin/python

import json
import os
import os.path
import re
import time

import matplotlib.pyplot as plt
import numpy as np
import poscar
from matplotlib.pyplot import MultipleLocator

import band
import pctl
import ppot
import scc_lib

# from pymatgen.ext.matproj import MPRester

charge_corr_dict = {'0Cs2AgInCl6': -0.10740636477538464,
                    '0Cs2HfCl6': -0.22865052200616062,
                    '1Cs2HfCl6': -0.22865052200616062 / 1.732,
                    '0Cs2KInCl6c2h': -0.06571622202097328,
                    '0Cs2KInCl6oh': -0.15055106860460557,
                    '1Cs2NaBiCl6': -0.1539912592779122,
                    '0Cs2NaBiCl6': -0.09500446458629402,
                    '0Cs2NaInCl6': -0.1968077259230447,
                    '0Cs2NaLaCl6': -0.15737025352611983,
                    '0Cs2NaScCl6': -0.1910407362901203,
                    '0Cs2NaYCl6': -0.1772608238927332,
                    '0Cs2SnCl6': -0.2562074178635485,
                    '1Cs2SnCl6': -0.2562074178635485 / 1.732,
                    '0Cs2ZrCl6': -0.22401592892112798,
                    '1Cs2ZrCl6': -0.22401592892112798 / 1.732,
                    '0KCl': -0.36247277151921625 / 2,
                    '0NaCl': -0.4186868438803664 / 2,
                    '0Rb3InCl6c2h': -0.08056157539366501,
                    '0Rb3InCl6oh': -0.13942761574705112,
                    '0Cs3Bi2Cl9': -0.1341969361103321,
                    "bao": -0.07760289267922081,
                    "cao": -0.05855249099089455,
                    "sao": -0.04836981934184487,
                    "cgo": -0.07150010428722961}

weight_for_each_queue = {"ckduan": 2, "smallopa": 1}


def recommend_queue(weight_for_each_queue, threhold):
    _, queue_list = jobid_in_queue()
    number_of_job = len(queue_list)
    while number_of_job > threhold:
        print(f"there are more than {threhold} jobs running, please wait ...")
        time.sleep(60)
        _, queue_list = jobid_in_queue()
        number_of_job = len(queue_list)
    recommend_dict = {queue: ((queue_list.count(queue) / number_of_job) if number_of_job > 0 else 0)
                             * weight_for_each_queue[queue]
                      for queue in weight_for_each_queue}

    argmin = np.array(list(recommend_dict.values())).argmin()
    return list(recommend_dict.keys())[argmin]


# 不定分组匹配
def match_times(sep, char, string):
    # char="\w" for element "\d" for nelement
    # 保持格式一致
    if sep == "\s":
        string = " " + string
    else:
        string = sep + string
    match_result = re.search("(" + sep + "+(" + char + "+))", string)
    match_list = []
    while match_result:
        unit = match_result.group(2)
        if char == "\d":
            unit = int(unit)
        match_list.append(unit)
        match_result = re.search("(" + sep + "+(" + char + "+)){" + str(len(match_list) + 1) + "}", string)
    return match_list


def jobid_in_queue():
    bjobs_return = scc_lib.linux_command("bjobs")[1:]
    # there are jobs
    if bjobs_return:
        map_func = lambda x: list(filter(None, x.strip().split()))[:4:3]
        id_list, queue_list = np.array(list(map(map_func, bjobs_return))).T
        return id_list.tolist(), queue_list.tolist()
    # there isn't jobs
    else:
        return [], []


def csv_split(string):
    return string.replace("\n", "").split(",,", 1)[0].split(",")


def origin_name(distort_name):
    origin_name = False
    distort_name = distort_name.replace("_half", "")
    if "p1" in distort_name:
        origin_name = distort_name.replace("p1", "").split("_attempt")[0]
    elif "m1" in distort_name:
        origin_name = distort_name.replace("m1", "").split("_attempt")[0]
    elif "ex" in distort_name:
        if "ex_ez" in distort_name:
            ex_name = "ex_ez"
        elif "ex_cz" in distort_name:
            ex_name = "ex_cz"
        elif "ex_cb" in distort_name:
            ex_name = "ex_cb"
        elif "ex_vb" in distort_name:
            ex_name = "ex_vb"
        elif "ex_nosym" in distort_name:
            ex_name = "ex_nosym"
        else:
            ex_name = "ex"
        origin_name = distort_name.replace(ex_name, "grd").split("_attempt")[0]
    return origin_name


def main_comp(spd_dict):
    max_value = 0
    max_key = False
    for i in spd_dict:
        if spd_dict[i] > max_value:
            max_value = spd_dict[i]
            max_key = i
    return max_key


# 交互/非交互访问material project
def matproj(formula, prop, path):
    # prop_list=['band_gap','cif','density','diel','e_above_hull','elasticity',
    #            'elements','energy','energy_per_atom',
    #            'formation_energy_per_atom','full_formula','hubbards',
    #            'icsd_id','icsd_ids','is_compatible','is_hubbard',
    #            'material_id','nelements','nsites','oxide_type','piezo',
    #            'pretty_formula','spacegroup','tags','task_ids',
    #            'total_magnetization','unit_cell_formula','volume']
    prop_list = ['band_gap', 'density', 'energy_per_atom',
                 'formation_energy_per_atom', 'material_id', 'full_formula',
                 'symbol', 'point_group', 'total_magnetization']
    point_group_dict = {"1": "C1",
                        "4": "C4",
                        "2": "C2",
                        "2/m": "C2h",
                        "mm2": "C2v",
                        "3": "C3",
                        "-6": "C3h",
                        "-3": "C3i",
                        "3m": "C3v",
                        "4/m": "C4h",
                        "4mm": "C4v",
                        "6": "C6",
                        "6/m": "C6h",
                        "6mm": "C6v",
                        "-1": "Ci",
                        "m": "Cs",
                        "222": "D2",
                        "-42m": "D2d",
                        "mmm": "D2h",
                        "32": "D3",
                        "-3m": "D3d",
                        "-6m2": "D3h",
                        "422": "D4",
                        "4/mmm": "D4h",
                        "622": "D6",
                        "6/mmm": "D6h",
                        "432": "O",
                        "m-3m": "Oh",
                        "-4": "S4",
                        "23": "T",
                        "-43m": "Td",
                        "m-3": "Th"}
    from pymatgen.ext.matproj import MPRester
    with MPRester("VUswQqBEWe4VFBZD25") as m:
        if prop:
            m.get_structure_by_material_id(
                sorted(m.get_data(formula), key=lambda x: x[prop])[0]["material_id"],
                final=True, conventional_unit_cell=True).to("poscar", path + "/CONTCAR")
            return False
        else:
            for i in range(len(prop_list)):
                print(i, prop_list[i])
            key = input("proprty to show: ")
            if key == "all":
                key = prop_list
            else:
                key = [prop_list[i] for i in match_times("\s", "\d", key)]

            value = []
            for i in m.get_data(formula):
                value.append([])
                for j in key:
                    if j == 'symbol':
                        value[-1].append(i['spacegroup'][j])
                    elif j == 'point_group':
                        value[-1].append(point_group_dict[i['spacegroup'][j]])
                    else:
                        value[-1].append(i[j])
            value.sort(key=lambda x: x[key.index('formation_energy_per_atom')])

            key.insert(0, "# ")
            for i in key:
                print("%-10s" % i[:10], end=" ")
            print()
            for i in value:
                i.insert(0, value.index(i))
                for j in i:
                    if type(j) == float:
                        j = round(j, 5)
                    print("%-10s" % j, end=" ")
                print()
            mp_list = match_times("\s", "\d", input("structure to get: "))
            if len(mp_list) == 1:
                m.get_structure_by_material_id(value[mp_list[0]][key.index("material_id")],
                                               final=True, conventional_unit_cell=True).to("poscar", path + "/CONTCAR")
                return False
            else:
                bundle_list = []
                for i in mp_list:
                    scc_lib.linux_command("mkdir " + path + "/" + value[i][key.index("material_id")])
                    m.get_structure_by_material_id(value[i][key.index("material_id")],
                                                   final=True, conventional_unit_cell=True).to("poscar",
                                                                                               path + "/" + value[i][
                                                                                                   key.index(
                                                                                                       "material_id")] + "/CONTCAR")
                    bundle_list.append(path + "/" + value[i][key.index("material_id")])
                return bundle_list


def plot_cc(energy_list, path):
    # energy_list 包含的能量依次为grd,ex,ex_at_grd,grd_at_ex

    c1 = 1
    c2 = 2
    num_size = 12
    word_size = 12

    grd, ex, ex_at_grd, grd_at_ex = energy_list
    c3 = 3
    plt.plot((c1, c1), (grd, ex_at_grd), 'ro--')
    plt.plot((c1, c2), (ex_at_grd, ex), 'ro--')
    plt.plot((c2, c2), (ex, grd_at_ex), 'ro--')
    plt.plot((c2, c1), (grd_at_ex, grd), 'ro--')
    plt.plot((c1, c2), (grd, ex), 'ro--')

    plt.annotate("grd", xy=(c1, grd), xytext=(-35, -5), textcoords='offset points', fontsize=word_size)
    plt.annotate("ex@grd", xy=(c1, ex_at_grd), xytext=(-60, -5), textcoords='offset points', fontsize=word_size)
    plt.annotate("ex", xy=(c2, ex), xytext=(10, -5), textcoords='offset points', fontsize=word_size)
    plt.annotate("grd@ex", xy=(c2, grd_at_ex), xytext=(10, -5), textcoords='offset points', fontsize=word_size)

    # 标记电荷转变能
    plt.annotate(str(round(ex_at_grd - grd, 2)) + "eV", xy=(c1, (grd + ex_at_grd) / 2), xytext=(-55, -5),
                 textcoords='offset points', fontsize=num_size)
    plt.annotate(str(round(ex_at_grd - ex, 2)) + "eV", xy=((c1 + c2) / 2, (ex + ex_at_grd) / 2), xytext=(10, 20),
                 textcoords='offset points', fontsize=num_size,
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    plt.annotate(str(round(ex - grd_at_ex, 2)) + "eV", xy=(c2, (ex + grd_at_ex) / 2), xytext=(10, -5),
                 textcoords='offset points', fontsize=num_size)
    plt.annotate(str(round(grd_at_ex - grd, 2)) + "eV", xy=((c1 + c2) / 2, (grd + grd_at_ex) / 2), xytext=(10, -30),
                 textcoords='offset points', fontsize=num_size,
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    plt.annotate(str(round(ex - grd, 2)) + "eV", xy=((c1 + c2) / 2, (grd + ex) / 2), xytext=(-20, -5),
                 textcoords='offset points', fontsize=num_size)

    x = np.linspace(2 * c1 - c2, 3 * c2 - 2 * c1, 100)
    coef = np.polyfit((c1, c2, c2 + 1e-5), (ex_at_grd, ex, ex), 2)
    plt.plot(x, np.poly1d(coef)(x))

    x = np.linspace(3 * c1 - 2 * c2, 2 * c2 - c1, 10)
    coef = np.polyfit((c1, c1 + 1e-5, c2), (grd, grd, grd_at_ex), 2)
    plt.plot(x, np.poly1d(coef)(x))

    plt.xlim(c1 - 3 * (c2 - c1), c2 + 3 * (c2 - c1))
    plt.xticks([])
    plt.xlabel('configurations')
    plt.ylim(grd - .6 * (ex - grd), ex + .6 * (ex - grd))
    plt.ylabel('energy(eV)')
    plt.title("configuration coordinate curves of " + path)
    # plt.legend()
    plt.savefig(path + ".png", dpi=300)
    plt.close('all')
    return None


replace_dict = {"ex_at_grd_pbe0": "sp@grd_w/o_soc",
                "ex_at_grd_pbe0_soc": "sp@grd",
                "ex_cb_at_grd_pbe0_soc": "cb@grd",
                "ex_vb_at_grd_pbe0_soc": "vb@grd",
                "ex_ez_pbe0_soc_temp1": "ez_static",
                "grd_at_ex_cz_pbe0_soc": "grd@cz",
                "grd_at_ex_ez_pbe0": "grd@ez_w/o_soc",
                "grd_at_ex_ez_pbe0_soc": "grd@ez",
                "grd_at_ex_cb_pbe0_soc": "grd@cb",
                "grd_at_ex_vb_pbe0_soc": "grd@vb",
                "relax_ex_cb_pbe0_soc": "cb",
                "relax_ex_cz_pbe0_soc": "cz",
                "relax_ex_ez_pbe0": "ez_w/o_soc",
                "relax_ex_ez_pbe0_soc": "ez",
                "relax_ex_vb_pbe0_soc": "vb",
                "relax_grd_pbe0": "grd_w/o_soc",
                "relax_grd_pbe0_soc": "grd",
                "relax_grdm1_pbe0_soc": "1+_+_cbm",
                "relax_grdp1_pbe0_soc": "1-_-_vbm",
                "relax_ex_cz_half_pbe0_soc": "cz_half",
                "relax_ex_ez_half_pbe0_soc": "ez_half",
                "ex_vb_at_grd_m1_pbe0_soc": "1vb@grd",
                "ex_vc_at_grd_m1_pbe0_soc": "1vc@grd",
                "ex_at_grd_m1_pbe0_soc": "1sp@grd",
                "ex_cb_at_grd_m1_pbe0_soc": "1cb@grd"}


def plot_cc3(energy_dict, path, xmin, xmax, y1min, y1max, y2min, y2max, rat):
    # energy_list 包含的能量依次为grd,ex,ex_at_grd,grd_at_ex

    word_size = 14

    # line_dict={"relax_grd":(("relax","ex_at"),("relax","ex_cb_at"),("relax","ex_vb_at"),("relax_grd","ex_vb_at_grd_m1"),("relax_grd","ex_cb_at_grd_m1"),("relax_grd","ex_at_grd_m1"),
    #                         ("grd","ex_ez"),("grd","ex_cz"),("grd","ex_cb"),("grd","ex_vb"),
    #                         ("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #            "ex_at_grd":(("ex_at_grd","relax_ex_ez"),("ex_at_grd","relax_ex_cz")),
    #            "ex_cb_at_grd":(("ex_cb_at_grd","relax_ex_cb"),()),
    #            "ex_vb_at_grd":(("ex_vb_at_grd","relax_ex_vb"),()),
    #            "relax_ex":(("relax","grd_at"),()),
    #            "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1"),()),
    #            "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
    #            "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}

    # fit_dict={"relax_grd":(("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #           "relax_ex":(("relax_ex_ez","ex_at_grd"),("relax_ex_cz","ex_at_grd")),
    #           "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd"),("relax_ex_cb","ex_cb_at_grd_m1")),
    #           "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd"),("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
    #           "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}
    for i in ["relax_grd_pbe0", "relax_ex_ez_pbe0", "grd_at_ex_ez_pbe0", "ex_at_grd_pbe0", "ex_ez_pbe0_soc_temp1",
              "ex_vb_at_grd_m1_pbe0_soc", "ex_vc_at_grd_m1_pbe0_soc", "ex_at_grd_m1_pbe0_soc",
              "ex_cb_at_grd_m1_pbe0_soc"]:
        if i in energy_dict:
            del energy_dict[i]
    print(energy_dict)
    # line_dict={"relax_grd":(("relax","ex_at"),("relax","ex_cb_at"),("relax","ex_vb_at"),("relax_grd","ex_vb_at_grd_m1"),("relax_grd","ex_cb_at_grd_m1"),("relax_grd","ex_at_grd_m1"),
    #                         ("grd","ex_ez"),("grd","ex_cz"),("grd","ex_cb"),("grd","ex_vb"),
    #                         ("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #            "ex_at_grd":(("ex_at_grd","relax_ex_ez"),("ex_at_grd","relax_ex_cz")),
    #            "ex_cb_at_grd":(("ex_cb_at_grd","relax_ex_cb"),()),
    #            "ex_vb_at_grd":(("ex_vb_at_grd","relax_ex_vb"),()),
    #            "relax_ex":(("relax","grd_at"),()),
    #            "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1"),()),
    #            "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
    #            "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}

    # fit_dict={"relax_grd":(("relax_grd","grd_at_ex_ez"),("relax_grd","grd_at_ex_cz"),("relax_grd","grd_at_ex_cb"),("relax_grd","grd_at_ex_vb")),
    #           "relax_ex":(("relax_ex_ez","ex_at_grd"),("relax_ex_cz","ex_at_grd")),
    #           "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd"),("relax_ex_cb","ex_cb_at_grd_m1")),
    #           "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd"),("relax_ex_vb","ex_vb_at_grd_m1")),
    #           "relax_ex_cz_pbe0_soc":(("relax_ex_cz_pbe0_soc","ex_at_grd_pbe0"),())}

    line_dict = {"relax_grd": (
        ("relax", "ex_at", "green"), ("relax", "ex_cb_at", "green"), ("relax", "ex_vb_at", "green"),
        ("relax_grd", "ex_vb_at_grd_m1"), ("relax_grd", "ex_cb_at_grd_m1"), ("relax_grd", "ex_at_grd_m1"),
        ("grd", "ex_ez"), ("grd", "ex_cz"), ("grd", "ex_cb"), ("grd", "ex_vb"),
        ("relax_grd", "grd_at_ex_ez"), ("relax_grd", "grd_at_ex_cz"), ("relax_grd", "grd_at_ex_cb"),
        ("relax_grd", "grd_at_ex_vb")),
        "ex_at_grd": (("ex_at_grd", "relax_ex_ez"), ("ex_at_grd", "relax_ex_cz")),
        "ex_cb_at_grd": (("ex_cb_at_grd", "relax_ex_cb"), ()),
        "ex_vb_at_grd": (("ex_vb_at_grd", "relax_ex_vb"), ()),
        # "relax_ex_cb":(("relax_ex_cb","ex_cb_at_grd_m1","red"),()),
        # "relax_ex_vb":(("relax_ex_vb","ex_vb_at_grd_m1"),("relax_ex_vb","ex_at_grd_m1")),
        # "relax_ex":(("relax","grd_at","blue"),()),
        "relax_ex_cb": (("relax_ex_cb", "ex_cb_at_grd_m1"), ("relax", "grd_at", "red"), ()),
        "relax_ex_vb": (
            ("relax_ex_vb", "ex_vb_at_grd_m1"), ("relax_ex_vb", "ex_at_grd_m1"), ("relax", "grd_at", "purple")),
        "relax_ex_ez": (("relax", "grd_at", "blue"), ()),
        "relax_ex_cz": (("relax", "grd_at", "cyan"), ()),
        "relax_ex_cz_pbe0_soc": (("relax_ex_cz_pbe0_soc", "ex_at_grd_pbe0"), ())}

    fit_dict = {"relax_grd": (("relax_grd", "grd_at_ex_ez", "green"), ("relax_grd", "grd_at_ex_cz", "green"),
                              ("relax_grd", "grd_at_ex_cb", "green"), ("relax_grd", "grd_at_ex_vb", "green")),
                "relax_ex": (("relax_ex_ez", "ex_at_grd", "blue"), ("relax_ex_cz", "ex_at_grd", "cyan")),
                "relax_ex_cb": (("relax_ex_cb", "ex_cb_at_grd", "red"), ("relax_ex_cb", "ex_cb_at_grd_m1")),
                "relax_ex_vb": (("relax_ex_vb", "ex_vb_at_grd", "purple"), ("relax_ex_vb", "ex_vb_at_grd_m1")),
                "relax_ex_cz_pbe0_soc": (("relax_ex_cz_pbe0_soc", "ex_at_grd_pbe0"), ())}
    ref = energy_dict["relax_grd_pbe0_soc"][1]

    f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

    for i in energy_dict:
        if type(energy_dict[i]) == tuple:
            energy_dict[i] = (energy_dict[i][0], energy_dict[i][1] - ref)
    for i in energy_dict:
        if type(energy_dict[i]) == tuple:
            for j in line_dict:
                if j in i:
                    for k in line_dict[j]:
                        if k:
                            if i.replace(*k[:2]) in energy_dict and type(energy_dict[i.replace(*k[:2])]) == tuple:
                                pos = list(zip(energy_dict[i], energy_dict[i.replace(*k[:2])]))
                                if abs(pos[0][0] - pos[0][1]) < 1e-3:
                                    # or abs(pos[1][0]-pos[1][1])<3
                                    ax.plot(*pos, linestyle="dashed", marker="o", color=k[2] if len(k) == 3 else "red")
                                    ax2.plot(*pos, linestyle="dashed", marker="o", color=k[2] if len(k) == 3 else "red")
                                # else:
                                #     x=np.linspace(*pos[0],2)
                                #     print([x,-(x-pos[0][0])/(pos[0][0]-pos[0][1])*(pos[1][1]-pos[1][0]-2.92)])
                                #     ax2.plot(x,-(x-pos[0][0])/(pos[0][0]-pos[0][1])*(pos[1][1]-pos[1][0]-2.92),'ro--')
                                #     print([x,-(x-pos[0][1])/(pos[0][0]-pos[0][1])*(pos[1][1]-pos[1][0]-2.92)+pos[1][1]])
                                #     ax.plot(x,-(x-pos[0][1])/(pos[0][0]-pos[0][1])*(pos[1][1]-pos[1][0]-2.92)+pos[1][1],'ro--')
                                if y1max < (pos[1][0] + pos[1][1]) / 2 < y2min:
                                    temp = str(round(abs(energy_dict[i.replace(*k[:2])][1] - energy_dict[i][1]), 2))
                                    temp = temp if len(temp) == 4 else temp + "0"

                                    if False:
                                        ax.annotate(temp,
                                                    xy=((pos[0][0] + pos[0][1]) / 2,
                                                        ((pos[1][0] + pos[1][1]) / 2 - 2) / 8 + y2min), xytext=(0, 0),
                                                    textcoords='offset points', fontsize=word_size)
                                        ax2.annotate(temp,
                                                     xy=((pos[0][0] + pos[0][1]) / 2,
                                                         ((pos[1][0] + pos[1][1]) / 2 - 2) / 8 + y2min), xytext=(0, 0),
                                                     textcoords='offset points', fontsize=word_size)
                                    else:
                                        ax.annotate(temp,
                                                    xy=((pos[0][0] + pos[0][1]) / 2,
                                                        ((pos[1][0] + pos[1][1]) / 2 - 2) / 8 + y1max), xytext=(0, 0),
                                                    textcoords='offset points', fontsize=word_size)
                                        ax2.annotate(temp,
                                                     xy=((pos[0][0] + pos[0][1]) / 2,
                                                         ((pos[1][0] + pos[1][1]) / 2 - 2) / 8 + y1max), xytext=(0, 0),
                                                     textcoords='offset points', fontsize=word_size)

                                # else:
                                #     ax.annotate(str(round(abs(energy_dict[i.replace(*k[:2])][1]-energy_dict[i][1]),2)),
                                #                  xy=((pos[0][0]+pos[0][1])/2,(pos[1][0]+pos[1][1])/2),xytext=(0,0),
                                #                  textcoords='offset points',fontsize=word_size)
                                #     ax2.annotate(str(round(abs(energy_dict[i.replace(*k[:2])][1]-energy_dict[i][1]),2)),
                                #                  xy=((pos[0][0]+pos[0][1])/2,(pos[1][0]+pos[1][1])/2),xytext=(0,0),
                                #                  textcoords='offset points',fontsize=word_size)
            for j in fit_dict:
                if j in i:
                    for k in fit_dict[j]:
                        if k:
                            if i.replace(*k[:2]) in energy_dict and type(energy_dict[i.replace(*k[:2])]) == tuple:
                                pos = list(zip(energy_dict[i], energy_dict[i.replace(*k[:2])]))
                                if False:
                                    pass
                                # if "relax_grd" in i:
                                #     if pos[0][1]>pos[0][0]:
                                #         x=np.linspace(pos[0][0]+1.5*abs(pos[0][0]-pos[0][1]),
                                #                       pos[0][0]
                                #                       ,50)
                                #     else:
                                #         x=np.linspace(pos[0][0],
                                #                       pos[0][0]-1.5*abs(pos[0][0]-pos[0][1])
                                #                       ,50)
                                else:
                                    x = np.linspace(pos[0][0] - 1.5 * abs(pos[0][0] - pos[0][1]),
                                                    pos[0][0] + 1.5 * abs(pos[0][0] - pos[0][1])
                                                    , 100)
                                coef = np.polyfit((pos[0][0] - 1e-5, pos[0][0] + 1e-5, pos[0][1]),
                                                  (pos[1][0], pos[1][0], pos[1][1]), 2)
                                ax.plot(x, np.poly1d(coef)(x), color=k[2] if len(k) == 3 else "red")
                                ax2.plot(x, np.poly1d(coef)(x), color=k[2] if len(k) == 3 else "red")
            # ax.plot(*energy_dict[i],"rx")
            # ax2.plot(*energy_dict[i],"rx")
            if not "at" in i:
                ax.annotate("" if replace_dict[i] == "grd" else replace_dict[i], xy=energy_dict[i], xytext=(0, 10),
                            textcoords='offset points', fontsize=word_size)
                ax2.annotate("" if replace_dict[i] == "grd" else replace_dict[i], xy=energy_dict[i], xytext=(0, 10),
                             textcoords='offset points', fontsize=word_size)
    # ax3=plt.gca()
    # ax.xaxis.set_major_locator(MultipleLocator(0.04))
    ax2.set_ylim(y1min, y1max)
    ax.set_ylim(y2min, y2max)
    ax.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop=False)
    ax2.xaxis.tick_bottom()
    d = .015
    plt.xlim(xmin, xmax)
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)

    # ax2.annotate("Bi",xy=(2.78,4),xytext=(0,0),textcoords='offset points',fontsize=14)
    # ax.annotate("Bi",xy=(2.79,4),xytext=(0,0),textcoords='offset points',fontsize=14)
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    # plt.subplots_adjust(left=0.125,bottom=0.1,right=0.9,top=0.9,wspace=0.2,hspace=0.35)
    plt.subplots_adjust(hspace=0.1)

    # plt.xlabel('Average bond length ($\AA$)',fontsize=16)
    plt.xlabel('Distortion ($\AA$)', fontsize=16)
    plt.ylabel('Energy (eV)', fontsize=16)
    ax.tick_params(labelsize=14)
    ax2.tick_params(labelsize=14)
    ax.set_aspect(rat)
    ax2.set_aspect(rat)
    plt.tight_layout()
    plt.savefig(path + ".png", dpi=300)
    plt.close('all')
    return None


# def plot_ctl(band_range, charge_list, path):
#     # charge transition level
#     x = np.linspace(band_range[0], band_range[1] + 0.5, 500)
#     # 不能给出交点
#     # plt.plot(x,[min([i[0]*j+i[1] for i in charge_list]) for j in x],color="red")
#
#     y = []
#     for i in range(len(charge_list)):
#         y.append(charge_list[i][0] * x + charge_list[i][1])
#
#     node_list = []
#     y_min = []
#     prev_min_key = False
#     for i in range(len(x)):
#         min_key = False
#         min_value = 1e5
#         for j in range(len(y)):
#             if y[j][i] < min_value:
#                 min_value = y[j][i]
#                 min_key = j
#         if type(prev_min_key) == bool:
#             prev_min_key = min_key
#         elif not min_key == prev_min_key:
#             node_list.append((x[i], (y[min_key][i] + y[prev_min_key][i]) * 0.5))
#             prev_min_key = min_key
#         y_min.append(min_value)
#
#     factor = 2e-4
#
#     y_min_max = max(y_min)
#     y_min_min = min(y_min)
#
#     plt.plot(x, y_min, color="red")
#     plt.plot([band_range[0], band_range[1]], [y_min[0], y_min[-1]], "ro")
#     plt.annotate(str(round(band_range[0], 2)) + "eV", xy=(band_range[0], y_min[0]), xytext=(-25, -15),
#                  textcoords='offset points', fontsize=14)
#     plt.annotate(str(round(band_range[1], 2)) + "eV", xy=(band_range[1], y_min[-1]), xytext=(-20, 15),
#                  textcoords='offset points', fontsize=14)
#     plt.plot([band_range[0], band_range[0]], [y_min_min * (1 + factor), y_min_max * (1 - factor)], linestyle='dashed',
#              color="blue")
#     plt.plot([band_range[1], band_range[1]], [y_min_min * (1 + factor), y_min_max * (1 - factor)], linestyle='dashed',
#              color="blue")
#
#     for i in node_list:
#         plt.plot(i[0], i[1], "ro", color="red")
#         plt.annotate(str(round(i[0], 2)) + "eV", xy=(i[0], i[1]), xytext=(-20, 10), textcoords='offset points',
#                      fontsize=14)
#         plt.plot([i[0], i[0]], [y_min_min * (1 + factor), y_min_max * (1 - factor)], linestyle='dashed', color="blue")
#
#     plt.xlabel('band(eV)')
#     plt.ylim(y_min_min * (1 + factor), )
#     plt.ylabel('energy(eV)')
#     plt.title("charge transition level of " + path)
#     # plt.legend()
#     plt.savefig(path + ".png", dpi=300)
#     plt.close('all')
#     return list(zip(x, y))


def plot_ctl3(split_list, title, path):
    interval = 1.5
    length = 1

    # 价带顶对齐
    # for i in split_list:
    #     avg=i["VBM"]
    #     for j in i:
    #         if not j=="system":
    #             i[j]-=avg

    # 真空能级对齐
    avg_dict = {"Cs2NaLaCl6": -5.584970968,
                "Cs2SnCl6": -6.227959018,
                "Cs2ZrCl6": -5.957924466,
                "Cs2NaInCl6": -6.0492546,
                "Cs2NaBiCl6": -5.844214162,
                "KCl": -4.221406774,
                "NaCl": -5.53,
                "Cs2AgInCl6": -6.633276321,
                "Cs2NaScCl6": -5.755407712,
                "Cs2NaYCl6": -5.407844399,
                "Rb3InCl6oh": -4.887988169,
                "Rb3InCl6c2h": -4.887988169,
                "Cs2HfCl6": -5.804700327,
                "Cs2KInCl6oh": -5.515516321,
                "Cs2KInCl6c2h": -5.515516321}

    for i in split_list:
        for j in i:
            if not j == "system":
                i[j] += avg_dict[i["system"][1:].split("_")[0]]

    line_list = []
    for i in split_list:
        temp_line_list = []
        for j in i:
            if not j == "system":
                temp_line_list.append([i[j], [], j, [0, 0]])
                # if j=="VBM":
                #     temp_line_list.append([i[j],[],i["system"],[0,-30]])
        line_list.append(temp_line_list)

    # line_list.sort(key=lambda x:max(list(zip(*x))[0]))
    for i in range(len(line_list)):
        for j in range(len(line_list[i])):
            if "VBM" in line_list[i][j][2]:
                plt.plot([i * interval, i * interval + length], [line_list[i][j][0], line_list[i][j][0]], color="black",
                         linewidth=1)
                plt.fill_between(np.linspace(i * interval - 0.25, i * interval + length + 0.25, 15), line_list[i][j][0],
                                 -10, color="black", alpha=0.5)
            elif "CBM" in line_list[i][j][2]:
                plt.plot([i * interval, i * interval + length], [line_list[i][j][0], line_list[i][j][0]], color="black",
                         linewidth=1)
                plt.fill_between(np.linspace(i * interval - 0.25, i * interval + length + 0.25, 15), line_list[i][j][0],
                                 2, color="black", alpha=0.25)
            elif line_list[i][j][2] in ["0/1", "-1/0"]:
                if line_list[i][j][0] > 1e4:
                    line_list[i][j][0] -= 1e5
                    plt.plot([i * interval, i * interval + length], [line_list[i][j][0], line_list[i][j][0]],
                             color="red", linestyle="--", linewidth=1)
                else:
                    plt.plot([i * interval, i * interval + length], [line_list[i][j][0], line_list[i][j][0]],
                             color="red", linewidth=1)
            elif "Ag" in line_list[i][j][2]:
                plt.plot([i * interval, i * interval + length], [line_list[i][j][0], line_list[i][j][0]], color="blue",
                         linewidth=1)
            else:
                plt.plot([i * interval, i * interval + length], [line_list[i][j][0], line_list[i][j][0]], color="cyan",
                         linewidth=1)
                # if line_list[i][j][2]:
                #     if "_" in line_list[i][j][2] and not "/" in line_list[i][j][2]:
                #         plt.annotate(line_list[i][j][2],xy=(i*interval,line_list[i][j][0]),
                #                   xytext=line_list[i][j][3],textcoords='offset points',fontsize=4,rotation="330")
                #     else:
                #         plt.annotate(line_list[i][j][2],xy=(i*interval,line_list[i][j][0]),
                #                   xytext=line_list[i][j][3],textcoords='offset points',fontsize=4,rotation="0")
    plt.xticks(np.linspace(0, len(split_list), len(split_list) + 1)[:-1] * interval + 0.5,
               [i["system"][1:].replace("Cs2", "").replace("Cl6", "").replace("_Bi", "").replace("_Sb", "") for i in
                split_list])
    plt.ylabel('energy (eV)')
    plt.xlim(-0.25, 14.75)
    plt.ylim(-9, 1)
    # plt.title("ctl")
    if path:
        plt.savefig(path + ".png", dpi=300)
    else:
        plt.show()
    plt.close("all")
    return None


def plot_nei(nll, path):
    fig, ax = plt.subplots(len(nll), sharex=True, sharey=True)
    for i in range(len(nll)):
        if i == 0:
            ax[i].set_title("neighbours of " + path)
        ax[i].set_yticks([])
        ax[i].annotate(nll[i][-1], xy=(3, 6), xytext=(80, -10), textcoords='offset points', fontsize=12)
        for j in nll[i][:-1]:
            ax[i].plot([j[1], j[1]], [0, j[2]], "r" if j[0] == "Cl" else "b", label=j[0])

    plt.subplots_adjust(hspace=0)
    plt.xlabel('distance(Å)')
    plt.savefig(path + ".png", dpi=300)
    plt.close("all")
    return None


def plot_pot(x, y, path):
    plt.plot(x[1:], y[1:], 'r')
    plt.xlabel(x[0])
    plt.ylabel(y[0])
    plt.title(path)
    plt.savefig(path + ".png", dpi=300)
    plt.close('all')
    return None


def plot_band(band, label, path):
    for i in band[1:]:
        if max(i) > -10:
            plt.plot(band[0], i, "r")
    for i in label:
        plt.xticks(label[0], label[1])
        plt.plot([label[0], label[0]], [-10, 10], "b--")

    # plt.xlabel(x[0])
    # plt.ylabel(y[0])
    # plt.title(path)
    plt.savefig(path + ".png", dpi=300)
    plt.close('all')
    return None


def read_procar(path):
    band_pattern = re.compile("band\s+(\d+) # energy\s+([\d\.\-]+) # occ.\s+([\d\.]+)")
    p = poscar.poscar()
    p.read(path + "/CONTCAR")
    spd_mode = False
    spd_begin = False
    spd_dict = {}
    current_band = {}
    band_list = []
    procar_list = []

    with open(path + "/PROCAR", "r") as f:
        for i in f.readlines():
            # band行
            if "occ." in i:
                current_band["num"], current_band["energy"], current_band["occ"] = [float(j) for j in
                                                                                    band_pattern.search(i).groups()]
                if current_band["num"] == 1 and band_list:
                    procar_list.append(band_list)
                    band_list = []
            # ion行
            elif " s " in i:
                if not spd_mode:
                    spd_mode = match_times("\s", "[\w\-]", i)
                    spd_pattern = re.compile("\s+".join(["([\d\.]+)" for j in range(len(spd_mode))]))
                spd_begin = True
            elif spd_begin:
                match_result = spd_pattern.search(i)
                if match_result:
                    # tot>0.01 获取主要部分
                    if float(match_result.groups()[-1]) > 0.01:
                        for j, k in zip(spd_mode[1:-1], match_result.groups()[1:-1]):
                            # 投影>0.001 获取主要部分
                            if float(k) > 0.001:
                                component = (p.label[int(match_result.groups()[0]) - 1][0], j.replace("x", "d")[0])
                                if not component in spd_dict:
                                    # spd_dict.append(component)
                                    spd_dict[component] = float(k)
                                else:
                                    spd_dict[component] += float(k)
                else:
                    # 当前band结束
                    spd_begin = False
                    if spd_dict:
                        for j in spd_dict:
                            spd_dict[j] = round(spd_dict[j], 3)
                    current_band["spd"] = spd_dict
                    spd_dict = {}
                    band_list.append(current_band)
                    current_band = {}
        procar_list.append(band_list)
    return procar_list


def energy(path):
    with open(path + "/OSZICAR", "r") as f:
        # match_result=re.search("E0=\s+([\+\-\.\dE]+)",f.readlines()[-1])
        # if match_result:
        #     return float(match_result.group(1))
        print(path, end=" ")
        last2lines = f.readlines()[-2:]
        if "E0" in last2lines[-1]:
            print(float(last2lines[-1].split("E0= ")[1].split(" ")[0]))
            return float(last2lines[-1].split("E0= ")[1].split(" ")[0])
        elif "E+" in last2lines[-1]:
            print(float([i for i in last2lines[-1].split(" ") if len(i) > 1 and i[:2] == "-0"][0]))
            return float([i for i in last2lines[-1].split(" ") if len(i) > 1 and i[:2] == "-0"][0])
        elif "E+" in last2lines[-2]:
            print(float([i for i in last2lines[-2].split(" ") if len(i) > 1 and i[:2] == "-0"][0]))
            return float([i for i in last2lines[-2].split(" ") if len(i) > 1 and i[:2] == "-0"][0])
        else:
            return False


class job:
    def __init__(self):
        pass

    def read_job(self, string):
        self.state = "wait"
        self.bundle = []
        detail_list = csv_split(string)
        template_dict = {
            "vasp": ["type", "structure", "wavecar", "occ", "spd1", "spd2", "chg",
                     "kmesh", "relax", "soc", "hdft", "nele", "spin", "ldau", "queue", "core", "encut", "gga", "mix",
                     "path"],
            "cc": ["type", "structure", "vib1", "vib2", "ini", "fin", "step", "path"],
            "cc2": ["type", "structure", "vib1", "vib2", "ini", "fin", "step", "path"],
            "para": ["type", "structure", "parameter", "ini", "fin", "step", "path"],
            "distort": ["type", "structure", "mode", "path"],
            "distort_sphere": ["type", "structure", "distort_list", "path"],
            "sub": ["type", "structure", "origin_site", "new_site", "path"],
            "matproj": ["type", "formula", "prop", "path"],
            "atom": ["type", "unitcell", "element", "path"],
            "slab": ["type", "structure", "vac_length", "path"],
            "slab2": ["type", "structure", "vac_length", "path"],
            "expand": ["type", "structure", "tm", "path"],
            "plot_split": ["type", "exm1", "exm1_soc", "grd_soc", "spd", "path"],
            "plot_split_cc": ["type", "cc_path", "spd", "path"],
            "plot_split_cc2": ["type", "cc1_path", "cc2_path", "spd", "path"],
            "plot_ctl": ["type", "unitcell", "pattern", "path"],
            "plot_ctl3": ["type", "p1", "m1", "unitcell", "path"],
            "ctl4": ["type", "system", "level"],
            "plot_nei": ["type", "atom", "range", "path"],
            "plot_cc": ["type", "grd", "ex", "ex_at_grd", "grd_at_ex", "path"],
            "plot_cc2": ["type", "mode", "path"],
            "plot_pot": ["type", "pot", "path"],
            "plot_band": ["type", "structure", "path"],
            "show_procar": ["type", "structure", "origin_structure", "spd1", "spd2"],
            "primcell": ["type", "structure", "path"],
            "formation": ["type", "formula", "defect_list", "host1", "host2", "name"]}
        self.detail = {}
        for i, j in zip(template_dict[detail_list[0]], detail_list):
            if j == "FALSE":
                self.detail[i] = False
            elif j == "TRUE":
                self.detail[i] = True
            elif i in ["ini", "fin", "spin", "encut", "vac_length"]:
                self.detail[i] = float(j)
            elif i in ["core", "range", "step", "relax"]:
                self.detail[i] = int(j)
            elif i == "nele":
                self.detail[i] = eval(j)
            elif i == "tm":
                self.detail[i] = []
                j = j.split("_")
                for k in range(len(j)):
                    if not k % 3:
                        self.detail[i].append([])
                    self.detail[i][-1].append(float(j[k]))
                self.detail[i] = np.array(self.detail[i])
            elif i == "ldau":
                self.detail[i] = {}
                j = j.split("_")
                for k in range(len(j)):
                    if not k % 4:
                        self.detail[i][j[k]] = []
                        ok = k
                    else:
                        self.detail[i][j[ok]].append(j[k])
            else:
                self.detail[i] = j
        return None

    def copy(self):
        ajob = job()
        ajob.state = self.state
        ajob.bundle = self.bundle
        ajob.detail = {}
        for i in self.detail:
            ajob.detail[i] = self.detail[i]
        return ajob

    def mk_kpoints_and_potcar(self):

        if self.detail["kmesh"] == "auto":
            scc_lib.linux_command("(echo 102;echo 2;echo 0.05)|vaspkit")
        elif self.detail["kmesh"] == "gamma":
            # scc_lib.linux_command("(echo 102;echo 2;echo 0)|vaspkit")
            scc_lib.linux_command("echo 'gam\n0\nGamma\n1 1 1\n0 0 0' >> KPOINTS")
            scc_lib.linux_command("(echo 103)|vaspkit")
        elif self.detail["kmesh"] == "line":
            scc_lib.linux_command("vaspkit -task 303")
            scc_lib.linux_command("cp KPATH.in KPOINTS")
        num_pattern = re.compile("([\d\.]+)")
        match_list = []
        for i in scc_lib.linux_command("grep 'parameters from PSCTR are:' POTCAR -B1"):
            match_result = num_pattern.search(i)
            if match_result:
                match_list.append(float(match_result.group(1)))
        return match_list

    def get_encut(self):
        encut_list = []
        with open("POTCAR", "r") as f:
            for i in f.readlines():
                if "ENMAX" in i:
                    encut_list.append(float(i.split("=")[1].split(";")[0]))
        return max(encut_list)

    def mk_incar(self):

        cwd = os.getcwd()
        os.chdir(self.detail["path"])
        p = poscar.poscar()
        p.read("POSCAR")
        nel = self.detail["nele"] + np.dot(p.number, self.mk_kpoints_and_potcar())
        incar = {
            "NCORE": int(np.power(self.detail["core"], 0.5) / 2) * 2,
            "NELECT": nel,
            "LREAL": "Auto",
            "ENCUT": self.get_encut() if self.detail["encut"] < 100 else self.detail["encut"],
            "ALGO": "All",
            "GGA": self.detail["gga"],
            "LORBIT": 11,
            "NSW": 0,
            "IBRION": -1,
            "EDIFF": 1e-4,
            "NELM": 50,
            "NELMIN": 5,
            "ISMEAR": 0,
            "SIGMA": 0.05
        }
        if not self.detail["wavecar"]:
            incar["NELMDL"] = -10
        if self.detail["relax"]:
            if type(self.detail["relax"]) == bool:
                incar["NSW"] = 25
            else:
                incar["NSW"] = self.detail["relax"]
            # if "ez" in self.detail["path"] or "cz" in self.detail["path"] or "temp" in self.detail["path"]:
            #     incar["IBRION"]=1
            # else:
            #     incar["IBRION"]=2
            incar["IBRION"] = 2
            # 计算介电张量
            if incar["NSW"] == 1:
                if "diele" in self.detail["path"]:
                    incar["IBRION"] = 8
                    incar["LEPSILON"] = ".TRUE."
                elif "pho" in self.detail["path"]:
                    incar["IBRION"] = 6
                    incar["NFREE"] = 2
                    incar["POTIM"] = 0.015

                incar["NCORE"] = 1
            else:
                incar["POTIM"] = 0.5
            incar["ISIF"] = 2
            incar["NELM"] = 30
            incar["EDIFFG"] = -0.02
            incar["ALGO"] = "Normal"
        if self.detail["soc"]:
            incar["LSORBIT"] = "True"
            incar["SAXIS"] = "0 0 1"
            incar["MAGMOM"] = "1000*0.5"
        if self.detail["hdft"]:
            incar["ALGO"] = "All"
            incar["LHFCALC"] = ".TRUE."
            incar["PRECFOCK"] = "Fast"
            if "pbe0" in self.detail["hdft"]:
                incar["AEXX"] = float("0." + self.detail["hdft"].replace("pbe0", ""))
            else:
                incar["AEXX"] = 0.25
                if self.detail["hdft"] == "hse03":
                    incar["HFSCREEN"] = 0.3
                elif self.detail["hdft"] == "hse06":
                    incar["HFSCREEN"] = 0.2
        if self.detail["chg"]:
            incar["ICHARG"] = 11
        os.chdir(cwd)

        if "slab" in self.detail["path"] or "atom" in self.detail["path"]:
            incar["LVHAR"] = ".TRUE."

        if self.detail["occ"]:
            incar["ISMEAR"] = -2
            if self.detail["occ"] == "easy":
                incar["FERWE"] = str(int(nel) - 1) + "*1 0 1 1000*0"
            else:
                fer = self.get_occ(0.5, 0.5) if "half" in self.detail["path"] else self.get_occ()
                if not type(fer) == bool:
                    if len(fer) > 1:
                        incar["FERWE"], incar["FERDO"] = fer
                    else:
                        incar["FERWE"] = fer[0]
                else:
                    if not fer:
                        incar["FERWE"] = "failed"

        os.chdir(self.detail["path"])
        if self.detail["mix"]:
            if self.detail["soc"]:
                incar["AMIX"] = 0.2
                incar["BMIX"] = 1e-5
                incar["AMIX_MAG"] = 0.8
                incar["BMIX_MAG"] = 1e-5
            else:
                incar["AMIX"] = 0.2
                incar["BMIX"] = 1e-5
        if not type(self.detail["spin"]) == bool:
            incar["ISPIN"] = 2
            incar["NUPDOWN"] = self.detail["spin"]
        elif self.detail["spin"]:
            incar["ISPIN"] = 2
        if self.detail["ldau"]:
            incar["LDAU"] = ".TRUE."
            incar["LDAUL"] = ""
            incar["LDAUU"] = ""
            incar["LDAUJ"] = ""
            for i in p.element:
                if i in self.detail["ldau"]:
                    for j, k in zip(self.detail["ldau"][i], ["LDAUL", "LDAUU", "LDAUJ"]):
                        incar[k] += j + " "
                else:
                    for j, k in zip(["-1", "0", "0"], ["LDAUL", "LDAUU", "LDAUJ"]):
                        incar[k] += j + " "
        with open("INCAR", "w") as f:
            for i in incar:
                f.write(i + "=" + str(incar[i]) + "\n")
        os.chdir(cwd)
        return None

    def runvasp(self, tag="normal"):

        cwd = os.getcwd()

        if tag == "rerun":
            pass
        elif tag == "normal":
            scc_lib.linux_command("mkdir " + self.detail["path"])
            scc_lib.linux_command("cp " + self.detail["structure"] + "/CONTCAR " + self.detail["path"] + "/POSCAR")
            if self.detail["wavecar"]:
                scc_lib.linux_command("cp " + self.detail["wavecar"] + "/WAVECAR " + self.detail["path"] + "/WAVECAR")
            if self.detail["chg"]:
                scc_lib.linux_command("cp " + self.detail["chg"] + "/CHGCAR " + self.detail["path"] + "/CHGCAR")

            self.mk_incar()
        os.chdir(self.detail["path"])
        hostname = scc_lib.linux_command("hostname")[0].replace("\n", "")
        if self.detail["soc"]:
            version = "soc"
        elif self.detail["kmesh"] == "gamma":
            version = "gam"
        else:
            version = "std"

        queue = recommend_queue(weight_for_each_queue, 75)
        print(f"0run -q {queue} -v {version}")
        scc_lib.linux_command(f"0run -q {queue} -v {version}")
        os.chdir(cwd)

        return None


class jobs:
    def __init__(self):

        self.job_list = []
        self.job_path_list = []
        for i in os.listdir():
            if ".csv" in i:
                with open(i, "r", encoding='UTF-8-sig') as f:
                    for j in f.readlines():
                        if not "#" in j:
                            ajob = job()
                            ajob.read_job(j)
                            # 去掉重复作业,不同策略除外
                            if not ajob.detail["type"] == "vasp":
                                if "path" in ajob.detail and ajob.detail["path"] in self.job_path_list:
                                    ajob = False
                            else:
                                for k in self.job_list:
                                    if k.detail == ajob.detail:
                                        ajob = False
                                        break
                            if ajob:
                                if "path" in ajob.detail:
                                    # print(ajob.detail["path"])
                                    self.job_path_list.append(ajob.detail["path"])
                                self.job_list.append(ajob)

    def find_job_by_path(self, path):
        i = False
        for i in self.job_list:
            if "path" in i.detail:
                if i.detail["path"] == path:
                    return i

    def update_vasp(self, path):

        global state
        ajob = self.find_job_by_path(path)
        listdir = os.listdir(path)
        running_list, _ = jobid_in_queue()
        # 获取最新的jobid
        if "jobid" in os.listdir(path):
            with open(path + "/jobid", "r") as f:
                jobids = f.readlines()

            jobid_list = []
            for line in jobids:
                if "<" in line:
                    jobid_list.append(line.split("<")[1].split(">")[0])
            if jobid_list:
                hostname = scc_lib.linux_command("hostname")[0].replace("\n", "")

                for jobid in jobid_list[::-1]:
                    if jobid in running_list:
                        state = "run"
                        print(f"{jobid} running")
                        break

                # finished
                else:
                    for jobid in jobid_list[::-1]:
                        # the last one
                        if jobid + ".log" in listdir:
                            print(f"{jobid} log")
                            with open(path + "/" + jobid + ".log", "r") as f:
                                log_lines = f.readlines()

                            if "writing wavefunctions" in log_lines[-1]:
                                if (ajob.detail["relax"] and not type(ajob.detail["relax"]) == int):
                                    if "reached required accuracy - stopping structural energy minimisation" in \
                                            log_lines[
                                                -2]:
                                        state = "done"
                                    else:
                                        print(f"rerun {path}")
                                        scc_lib.linux_command("cp " + path + "/CONTCAR " + path + "/POSCAR")
                                        ajob.runvasp("rerun")
                                        state = "run"
                                else:
                                    state = "done"
                            else:
                                state = "fail"
                            break
                    else:
                        state = "fail"
                        print(f"remove {path}")
                        # print(to_be_removed)
            else:
                print("empty jobid")
                state = "fail"

        else:
            print(f"resub {path}")
            state = "fail"
        # print(f"{state} {path}")
        return state

    def update(self, file, delay):

        jobid_template = re.compile("(\d+)")
        for i in self.job_list:
            # print(i.detail)
            # print(i.bundle)
            # print(i.state)
            if i.detail["type"] == "vasp":

                # 判断vasp计算中的POSCAR对称性
                # if os.path.isdir(i.detail["path"]):
                #     pwd=os.getcwd()
                #     os.chdir(i.detail["path"])
                #     for j in scc_lib.linux_command("vaspkit -task 601"):
                #         if "Point Group" in j:
                #             file.write("point "+j.split("[")[1].split("]")[0]+" "+i.detail["path"]+"\n")
                #             break
                #     os.chdir(pwd)

                if i.state in ["wait", "run"]:
                    # 更新vasp的bundle
                    job_structure = self.find_job_by_path(i.detail["structure"])
                    if i.detail["path"] in self.job_path_list:
                        # 不是子任务
                        if job_structure.bundle:
                            # bundle
                            if i.bundle:
                                state_list = [self.find_job_by_path(j).state for j in i.bundle]
                                if "fail" in state_list:
                                    i.state = "fail"
                                elif set(state_list) == {"done"}:
                                    i.state = "done"
                            elif os.path.isdir(i.detail["path"]):
                                # bundle丢失
                                for j in os.listdir(i.detail["structure"]):
                                    self.job_list.append(i.copy())
                                    self.job_list[-1].detail["structure"] = self.job_list[-1].detail[
                                                                                "structure"] + "/" + j
                                    self.job_list[-1].detail["path"] = self.job_list[-1].detail["path"] + "/" + j
                                    if self.job_list[-1].detail["wavecar"]:
                                        self.job_list[-1].detail["wavecar"] = self.job_list[-1].detail[
                                                                                  "wavecar"] + "/" + j

                                    # self.job_list[-1].state=self.update_vasp(self.job_list[-1].detail["path"])
                                    i.bundle.append(self.job_list[-1].detail["path"])
                                    i.state = "run"

                            elif (job_structure.state == "done"
                                  and (not i.detail["wavecar"] or self.find_job_by_path(
                                        i.detail["wavecar"]).state == "done")
                                  and (type(i.detail["occ"]) == bool or i.detail[
                                        "occ"] == "easy" or self.find_job_by_path(
                                        i.detail["occ"]).state == "done")
                                  and (not i.detail["chg"] or self.find_job_by_path(i.detail["chg"]).state == "done")):
                                # bundle未提交
                                scc_lib.linux_command("mkdir " + i.detail["path"])
                                for j in os.listdir(i.detail["structure"]):
                                    self.job_list.append(i.copy())
                                    self.job_list[-1].detail["structure"] = self.job_list[-1].detail[
                                                                                "structure"] + "/" + j
                                    self.job_list[-1].detail["path"] = self.job_list[-1].detail["path"] + "/" + j
                                    if self.job_list[-1].detail["wavecar"]:
                                        if os.path.isdir(self.job_list[-1].detail["wavecar"] + "/" + j):
                                            self.job_list[-1].detail["wavecar"] = self.job_list[-1].detail[
                                                                                      "wavecar"] + "/" + j
                                    self.job_list[-1].runvasp()
                                    self.job_list[-1].state = "run"
                                    i.bundle.append(self.job_list[-1].detail["path"])
                                    i.state = "run"
                        else:
                            # no bundle
                            if os.path.isdir(i.detail["path"]):
                                i.state = self.update_vasp(i.detail["path"])
                            elif (job_structure.state == "done"
                                  and (not i.detail["wavecar"] or self.find_job_by_path(
                                        i.detail["wavecar"]).state == "done")
                                  and (type(i.detail["occ"]) == bool or i.detail[
                                        "occ"] == "easy" or self.find_job_by_path(i.detail["occ"]).state == "done")
                                  and (not i.detail["chg"] or self.find_job_by_path(i.detail["chg"]).state == "done")):
                                i.runvasp()
                                i.state = "run"
                    # 路径存在
                    elif os.path.isdir(i.detail["structure"]):
                        # 子任务
                        if os.path.isdir(i.detail["path"]):
                            i.state = self.update_vasp(i.detail["path"])
                        else:
                            i.runvasp()
                            i.state = "run"

            # 非vasp
            elif i.state == "wait":
                # 1.prepared?
                # 2.作业路径
                # 3.bundle?
                # 4.p.read
                # 5.done
                if i.detail["type"] == "matproj":
                    if os.path.isdir(i.detail["path"]):
                        if not "CONTCAR" in os.listdir(i.detail["path"]):
                            for j in os.listdir(i.detail["path"]):
                                i.bundle.append(i.detail["path"] + "/" + j)
                        i.state = "done"
                    else:
                        scc_lib.linux_command("mkdir " + i.detail["path"])
                        i.bundle = matproj(i.detail["formula"], i.detail["prop"], i.detail["path"])
                        i.state = "done"

                elif i.detail["type"] == "atom":
                    if os.path.isdir(i.detail["path"]):
                        i.state = "done"
                    else:
                        p = poscar.poscar()
                        p.read(f"{i.detail['unitcell']}/CONTCAR")
                        scc_lib.linux_command("mkdir " + i.detail["path"])
                        p.atom(i.detail["element"]).write(i.detail["path"] + "/CONTCAR")
                        i.state = "done"


                elif i.detail["type"] == "formation" and not delay:
                    cp = []
                    with open("js_" + i.detail["formula"], "r") as f:
                        for j in f.readlines():
                            cp.append(json.loads(j.replace("\n", "")))
                    for j in cp:
                        j["chempot"]["vac"] = 0

                    # pot = {"Cs": 0.388, "Na": 0.142, "Bi": 0.303, "Cl": 0.15, "vac": 0, "O": 0.059, "Sn": 0.294,
                    #        "Ag": 0.189, "In": 0.275, "Sb": 0.29}

                    # pot = {"Bi": 0.087,
                    #        "Cl": 0.043,
                    #        "Cs": 0.094,
                    #        "Sb": 0.079,
                    #        "Se": 0.048,
                    #        "Sn": 0.081,
                    #        "Te": 0.081,
                    #        "Zr": 0.083,
                    #        "Hf": 0.055,
                    #        "vac": 0}
                    # pot = {"Pb": 0.161,
                    #        "O": 0.032,
                    #        "Ca": 0.131,
                    #        "Bi": 0.169,
                    #        "Al": 0.109,
                    #        "Fe": 0.074,
                    #        "Mn": 0.091,
                    #        "Sr": 0.147,
                    #        "Ba": 0.175,
                    #        "vac": 0}
                    # pot={"Cs":0.372,"Na":0.15,"Bi":0.312,"Cl":0.151,"vac":0}
                    # new
                    # pot={"Cs": -0.084665179,"Na": -0.035525564,"Bi": -0.078202227,"Cl": -0.037649451,"vac":0}

                    # for j in pot:
                    #     pot[j]/=4
                    path_dict = json.load(open("path.json", "r"))
                    pot = {"Ca": 0.222, "Ga": 0.153, "Bi": 0.278, "O": 0.054, "vac": 0}
                    energy_dict = {}
                    ajob = self.find_job_by_path(path_dict["band"])
                    ajob_soc = self.find_job_by_path(path_dict["soc"])
                    procar = read_procar(ajob.detail["path"])[0]
                    for j in procar:
                        if j["occ"] < 0.3:
                            cbm1 = j["energy"]
                            vbm1 = procar[procar.index(j) - 1]["energy"]
                            break
                    else:
                        assert False, "band gap not clear."

                    grd_energy = energy(ajob.detail["path"]) - charge_corr_dict[ajob.detail["path"].split("/")[0]] * \
                                 ajob.detail["nele"] ** 2

                    grd_energy_soc = energy(ajob_soc.detail["path"]) - charge_corr_dict[
                        ajob_soc.detail["path"].split("/")[0]] * \
                                     ajob.detail["nele"] ** 2

                    for j_path in path_dict["ctl"]:
                        j = self.find_job_by_path(j_path)
                        print(j_path)
                        temp_job = j
                        while True:
                            asub = self.find_job_by_path(temp_job.detail["structure"])
                            if asub.detail["type"] == "sub":
                                break
                            else:
                                temp_job = asub

                        origin_site = asub.detail["origin_site"].split("-")
                        new_site = asub.detail["new_site"].split("-")
                        defect_name = j.detail["path"].split("_")[1] + "_" + j.detail["path"].split("_")[-1]
                        if not energy_dict or not defect_name in energy_dict:
                            energy_dict[defect_name] = []
                        temp_energy_dict = {"defect": energy(j.detail["path"]),
                                            "madebula": -charge_corr_dict[j.detail["path"].split("/")[0]],
                                            "qq": -j.detail["nele"],
                                            "host": grd_energy_soc if j.detail["soc"] else grd_energy, "potalign": 0}
                        temp_str1 = ""
                        for k in origin_site:
                            temp_energy_dict["potalign"] -= pot[k]
                            temp_str1 += "-" + k + f":{pot[k]:.3f}"

                        for k in new_site:
                            temp_energy_dict["potalign"] += pot[k]
                            temp_str1 += "+" + k + f":{pot[k]:.3f}"

                        temp_energy_dict["potalign"] *= temp_energy_dict["qq"]
                        temp_energy_dict["img"] = temp_energy_dict["madebula"] * temp_energy_dict["qq"] * \
                                                  temp_energy_dict["qq"]

                        temp_energy_dict["energy"] = temp_energy_dict["defect"] - temp_energy_dict["host"] + \
                                                     temp_energy_dict["img"] + temp_energy_dict["potalign"]

                        temp_str = "_energy_=_defect_-_host_+(_img_=_madebula_*_qq_^2)+(_potalign_=" + temp_str1 + "*_qq_)"
                        temp_str1 = ""
                        for k in temp_str.split("_"):
                            if k in temp_energy_dict:
                                temp_str1 += k + f":{temp_energy_dict[k]:.3f}"
                                pass
                            else:
                                temp_str1 += k

                        # temp_str=f'energy:{temp_energy_dict["energy"]:.3f}=defect:{temp_energy_dict["defect"]:.3f}-host:{temp_energy_dict["host"]:.3f}'
                        # temp_str+=f'+(img:{temp_energy_dict["img"]:.3f}=madebula:{temp_energy_dict["madebula"]:.3f}*qq:{temp_energy_dict["qq"]}^2)'
                        # temp_str+=f'+(potalign:{temp_energy_dict["potalign"]:.3f}=('+temp_str1+f')*qq:{temp_energy_dict["qq"]}^2)'
                        # 打印计算过程
                        print(temp_str1)
                        energy_dict[defect_name].append(
                            [-j.detail["nele"], temp_energy_dict["energy"], defect_name, origin_site, new_site])

                    fermi_where = {"BiCaO3": 0.6969405148484848,
                                   "BiGa": 2.4847444442424242,
                                   "Bi2CaO4": 1.5756915987878788}
                    print(energy_dict)
                    for k in cp:
                        # for defect in ["BiCa", "BiGa", "CaGa", "GaCa", "vCa", "vGa", "vO"]:
                        charged_energy_dict = {j: [
                            [l[0],
                             l[1] + (sum([k["chempot"][m] for m in l[3]]) - sum([k["chempot"][m] for m in l[4]]))]
                            # for l in energy_dict[j]] for j in energy_dict if defect in j}
                            for l in energy_dict[j]] for j in energy_dict}
                        print(i.detail["name"] + i.detail["formula"] + "_".join(k["tag"]))

                        print(charged_energy_dict)

                        pctl.plot_ctl(charged_energy_dict,
                                      150,
                                      i.detail["name"] + i.detail["formula"] + "_".join(k["tag"]),
                                      [vbm1, cbm1, 100], fermi_where[k["tag"][-1].split("_")[-1]])
                    i.state = "done"

                elif i.detail["type"] == "show_procar":
                    if (self.find_job_by_path(i.detail["structure"]).state == "done"
                            and self.find_job_by_path(i.detail["origin_structure"]).state == "done"):

                        job_structure = self.find_job_by_path(i.detail["structure"])
                        file.write(i.detail["structure"] + "\n")
                        # 和get_occ 相同
                        procar = read_procar(i.detail["structure"])
                        start = False
                        end = False
                        if i.detail["spd1"]:
                            if "_" in i.detail["spd1"]:
                                i.detail["spd1"] = tuple(i.detail["spd1"].split("_"))
                            for j in procar[-1]:
                                if j["occ"] <= 0.5:
                                    # if type(i.detail["spd1"])==tuple:
                                    #     if not i.detail["spd1"] in j["spd"]:
                                    #         job_structure.state="wrong"
                                    # else:
                                    #     if not j["spd"] or not i.detail["spd1"] in list(zip(*j["spd"].keys()))[0]:
                                    #         job_structure.state="wrong"
                                    # start=procar[-1].index(j)

                                    if not "no" in i.detail["spd1"]:
                                        if not i.detail["spd1"] in j["spd"]:
                                            job_structure.state = "wrong"
                                    else:
                                        if j["spd"]:
                                            for k in j["spd"]:
                                                if i.detail["spd1"][1] in k and j["spd"][k] > 0.05:
                                                    job_structure.state = "wrong"
                                                    break

                                    start = procar[-1].index(j)
                                    break

                        if i.detail["spd2"]:
                            if "_" in i.detail["spd2"]:
                                i.detail["spd2"] = tuple(i.detail["spd2"].split("_"))
                            for j in procar[0][::-1]:
                                if j["occ"] >= 0.5:
                                    # if type(i.detail["spd2"])==tuple:
                                    #     if not i.detail["spd2"] in j["spd"]:
                                    #         job_structure.state="wrong"
                                    # else:
                                    #     if not j["spd"] or not i.detail["spd2"] in list(zip(*j["spd"].keys()))[0]:
                                    #         job_structure.state="wrong"

                                    if not "no" in i.detail["spd2"]:
                                        if not i.detail["spd2"] in j["spd"]:
                                            job_structure.state = "wrong"
                                    else:
                                        if j["spd"]:
                                            for k in j["spd"]:
                                                if i.detail["spd2"][1] in k and j["spd"][k] > 0.05:
                                                    job_structure.state = "wrong"
                                                    break

                                    end = procar[0].index(j)
                                    break

                        real_start = start - 5 if start else end - 5
                        end = end + 5 if end else start + 5

                        # 检查结构
                        distort_poscar = poscar.poscar()
                        origin_poscar = poscar.poscar()
                        origin_path = i.detail["origin_structure"]

                        # 根据预测的畸变计算控制占据
                        if not "attempt" in i.detail["structure"] and "relax" in i.detail["structure"]:
                            if i.detail["structure"].replace("relax_", "") + "_attempt" in self.job_path_list:
                                if (job_structure.state in ["wrong", "fail"] or
                                        self.find_job_by_path(job_structure.detail["structure"]).state in ["wrong",
                                                                                                           "fail"]):
                                    ajob = self.find_job_by_path(
                                        i.detail["structure"].replace("relax_", "") + "_attempt")
                                    ajob.detail["structure"] = ajob.detail["structure"].replace("_attempt", "")
                                elif job_structure.state in ["done", "run", "wait"]:
                                    for j in self.job_path_list:
                                        if (i.detail["structure"].replace("relax_", "") + "_attempt" in j or
                                                i.detail["structure"] + "_attempt" in j):
                                            if os.path.isdir(j):
                                                if "jobid" in os.listdir(j):
                                                    with open(j + "/jobid", "r") as f:
                                                        scc_lib.linux_command(
                                                            "bkill " + jobid_template.search(f.readlines()[-1]).group(
                                                                1))
                                                scc_lib.linux_command("mv " + j + " " + j + "_nomore")
                                                self.find_job_by_path(j).state = "wait"

                        origin_poscar.read(origin_path + "/CONTCAR")
                        distort_poscar.read(i.detail["structure"] + "/CONTCAR")
                        distortion = distort_poscar - origin_poscar
                        # 对于低对称情况缺陷也可能移动,故不能直接使用poscar.distance()
                        nearest_distortion = np.linalg.norm(
                            distortion.position[
                                distort_poscar.label.index(
                                    distort_poscar.distance(distort_poscar.label[0])[1][0])].dot(distortion.axis))

                        distort_list = []
                        current_element = False
                        neis = 0
                        ele_copy = False
                        if "unitcell" in i.detail["structure"]:
                            ban_list = ["Br", "F", "Cl", "Cs", "Na", "K", "Ag"]
                            ele_copy = distort_poscar.element[:]
                            for j in ban_list:
                                if len(ele_copy) == 1:
                                    break
                                elif j in ele_copy:
                                    ele_copy.remove(j)
                        ele_copy = (ele_copy[0], 1) if ele_copy else distort_poscar.label[0]
                        for j in distort_poscar.distance(ele_copy):
                            if not current_element == j[0][0]:
                                neis += 1
                                if neis > 3:
                                    break
                                else:
                                    current_element = j[0][0]
                            temp_distort_list = []
                            file.write("%-12s" % str(j[0]))
                            temp_distort_list.append(j[0][0])
                            abs_origin = (origin_poscar.position[distort_poscar.label.index(j[0])]
                                          - origin_poscar.position[distort_poscar.label.index(ele_copy)])
                            for k in range(len(abs_origin)):
                                if abs_origin[k] > 0.5:
                                    abs_origin[k] -= 1
                                elif abs_origin[k] <= -0.5:
                                    abs_origin[k] += 1
                            abs_origin = abs_origin.dot(origin_poscar.axis)
                            for k in abs_origin:
                                file.write("%-8s" % str(round(k, 3)))

                            file.write("%-8s" % str(round(np.linalg.norm(abs_origin), 3)))
                            temp_distort_list.append(round(np.linalg.norm(abs_origin), 3))
                            abs_distortion = distortion.position[distort_poscar.label.index(j[0])].dot(distortion.axis)
                            for k in abs_distortion:
                                file.write("%-8s" % str(round(k / nearest_distortion, 1)))

                            file.write("%-8s" % str(round(np.linalg.norm(abs_distortion), 3)))

                            if np.linalg.norm(abs_origin):
                                file.write(
                                    "%-8s" % str(round(abs_origin.dot(abs_distortion) / np.linalg.norm(abs_origin), 3)))
                                temp_distort_list.append(
                                    round(abs_origin.dot(abs_distortion) / np.linalg.norm(abs_origin), 3))
                            else:
                                file.write("%-8s" % str(round(np.linalg.norm(abs_distortion), 3)))
                            # if not temp_distort_list in distort_list:
                            distort_list.append(temp_distort_list)
                            file.write("\n")
                        del distort_list[0]
                        file.write(job_structure.state + " distort " + i.detail["structure"] + " ")
                        if distort_list[0] == distort_list[5]:
                            distort_mode = "A"
                        elif distort_list[0] == distort_list[3] and distort_list[4] == distort_list[5]:
                            distort_mode = "Ee"
                        elif distort_list[0] == distort_list[1] and distort_list[2] == distort_list[5]:
                            distort_mode = "Ec"
                        else:
                            distort_mode = "unk"
                        file.write(distort_mode + " ")
                        prev = False
                        for j in distort_list:
                            if not prev == j:
                                for k in j:
                                    file.write(str(k) + " ")
                                prev = j
                        file.write("\n")

                        # attempt成功则替换
                        if ("attempt" in i.detail["structure"]
                                and "relax" in i.detail["structure"]
                                and job_structure.state == "done"
                                and self.find_job_by_path(i.detail["structure"].split("_attempt")[0]).state in ["fail",
                                                                                                                "wrong"]):
                            scc_lib.linux_command("mv " + i.detail["structure"].split("_attempt")[0]
                                                  + " " + i.detail["structure"].split("_attempt")[0] + "_fail")
                            scc_lib.linux_command(
                                "mv " + i.detail["structure"] + " " + i.detail["structure"].split("_attempt")[0])
                            self.find_job_by_path(i.detail["structure"].split("_attempt")[0]).state = "done"

                        # 检查激发态
                        for j in range(end)[real_start:]:
                            for k in procar:
                                file.write(str(k[j]) + "\n")

                        i.state = "done"

                    elif ("relax" in i.detail["structure"]
                          and not "attempt" in i.detail["structure"]
                          and self.find_job_by_path(i.detail["structure"]).state == "fail"
                          and i.detail["structure"].replace("relax_", "") + "_attempt" in self.job_path_list):
                        ajob = self.find_job_by_path(i.detail["structure"].replace("relax_", "") + "_attempt")
                        ajob.detail["structure"] = ajob.detail["structure"].replace("_attempt", "")

                        i.state = "done"


                elif i.detail["type"] == "sub":
                    if os.path.isdir(i.detail["path"]):
                        listdir = os.listdir(i.detail["path"])
                        if not "CONTCAR" in listdir:
                            for j in listdir:
                                i.bundle.append(i.detail["path"] + "/" + j)
                        i.state = "done"

                    elif self.find_job_by_path(i.detail["structure"]).state == "done":
                        p = poscar.poscar()
                        if self.find_job_by_path(i.detail["structure"]).bundle:
                            sub_list = []
                            for j in self.find_job_by_path(i.detail["structure"]).bundle:
                                p.read(j + "/CONTCAR")
                                sub_list += p.sub(i.detail["origin_site"], i.detail["new_site"])
                            if not os.path.isdir(i.detail["path"]):
                                scc_lib.linux_command("mkdir " + i.detail["path"])
                            if len(sub_list) == 1:
                                sub_list[0][1].write(i.detail["path"] + "/CONTCAR")
                            else:
                                file.write("\n" + i.detail["path"] + "\n")
                                for j in sub_list:
                                    file.write("sub " + str(j[0]) + "\n")
                                    for k in j[1].distance(j[1].label[0])[:20]:
                                        file.write(str(k) + "\n")
                                    if not os.path.isdir(i.detail["path"] + "/" + str(sub_list.index(j))):
                                        scc_lib.linux_command(
                                            "mkdir " + i.detail["path"] + "/" + str(sub_list.index(j)))
                                        j[1].write(i.detail["path"] + "/" + str(sub_list.index(j)) + "/CONTCAR")
                                    i.bundle.append(i.detail["path"] + "/" + str(sub_list.index(j)))
                            i.state = "done"

                        else:
                            p.read(i.detail["structure"] + "/CONTCAR")
                            if i.detail["origin_site"] == "vac":
                                sub_list = p.inter(i.detail["new_site"])
                            else:
                                sub_list = p.sub(i.detail["origin_site"], i.detail["new_site"])
                            if not os.path.isdir(i.detail["path"]):
                                scc_lib.linux_command("mkdir " + i.detail["path"])
                            if len(sub_list) == 1:
                                sub_list[0][1].write(i.detail["path"] + "/CONTCAR")
                            else:
                                file.write("\n" + i.detail["path"] + "\n")
                                for j in sub_list:
                                    file.write("sub " + str(j[0]) + "\n")
                                    for k in j[1].distance(j[1].label[0])[:20]:
                                        file.write(str(k) + "\n")
                                    if not os.path.isdir(i.detail["path"] + "/" + str(j[0])):
                                        scc_lib.linux_command("mkdir " + i.detail["path"] + "/" + str(j[0]))
                                        j[1].write(i.detail["path"] + "/" + str(j[0]) + "/CONTCAR")
                                    i.bundle.append(i.detail["path"] + "/" + str(j[0]))
                            i.state = "done"

                elif i.detail["type"] == "cc":
                    if os.path.isdir(i.detail["path"]):
                        for j in os.listdir(i.detail["path"]):
                            i.bundle.append(i.detail["path"] + "/" + j)
                        i.state = "done"
                    if (self.find_job_by_path(i.detail["structure"]).state == "done" and
                            self.find_job_by_path(i.detail["vib1"]).state == "done" and
                            self.find_job_by_path(i.detail["vib2"]).state == "done"):
                        p = poscar.poscar()
                        v1 = poscar.poscar()
                        v2 = poscar.poscar()
                        p.read(i.detail["structure"] + "/CONTCAR")
                        v1.read(i.detail["vib1"] + "/CONTCAR")
                        v2.read(i.detail["vib2"] + "/CONTCAR")
                        if not os.path.isdir(i.detail["path"]):
                            scc_lib.linux_command("mkdir " + i.detail["path"])
                        for j in p.cc(v1 - v2, i.detail["ini"], i.detail["fin"], i.detail["step"]):
                            if not os.path.isdir(i.detail["path"] + "/cc_" + str(j[0])):
                                scc_lib.linux_command("mkdir " + i.detail["path"] + "/cc_" + str(j[0]))
                                j[1].write(i.detail["path"] + "/cc_" + str(j[0]) + "/CONTCAR")
                            i.bundle.append(i.detail["path"] + "/cc_" + str(j[0]))
                        i.state = "done"

                elif i.detail["type"] == "cc2":
                    if os.path.isdir(i.detail["path"]):
                        for j in os.listdir(i.detail["path"]):
                            i.bundle.append(i.detail["path"] + "/" + j)
                        i.state = "done"
                    if (self.find_job_by_path(i.detail["structure"]).state == "done" and
                            self.find_job_by_path(i.detail["vib1"]).state == "done" and
                            self.find_job_by_path(i.detail["vib2"]).state == "done"):
                        p = poscar.poscar()
                        v1 = poscar.poscar()
                        v2 = poscar.poscar()
                        p.read(i.detail["structure"] + "/CONTCAR")
                        v1.read(i.detail["vib1"] + "/CONTCAR")
                        v2.read(i.detail["vib2"] + "/CONTCAR")
                        if not os.path.isdir(i.detail["path"]):
                            scc_lib.linux_command("mkdir " + i.detail["path"])
                        for j in p.cc(v1 - p, i.detail["ini"], i.detail["fin"], i.detail["step"]):
                            for k in j[1].cc(v2 - p, i.detail["ini"], i.detail["fin"], i.detail["step"]):
                                if not os.path.isdir(i.detail["path"] + "/cc_" + str(j[0]) + "_" + str(k[0])):
                                    scc_lib.linux_command(
                                        "mkdir " + i.detail["path"] + "/cc_" + str(j[0]) + "_" + str(k[0]))
                                    k[1].write(i.detail["path"] + "/cc_" + str(j[0]) + "_" + str(k[0]) + "/CONTCAR")
                                i.bundle.append(i.detail["path"] + "/cc_" + str(j[0]) + "_" + str(k[0]))
                        # for j in p.cc(v1-v2,i.detail["ini"],i.detail["fin"],i.detail["step"]):
                        #     for k in j[1].cc(p.cc2(),0,0.5,6):
                        #         if not os.path.isdir(i.detail["path"]+"/cc_"+str(j[0])+"_"+str(k[0])):
                        #             scc_lib.linux_command("mkdir "+i.detail["path"]+"/cc_"+str(j[0])+"_"+str(k[0]))
                        #             k[1].write(i.detail["path"]+"/cc_"+str(j[0])+"_"+str(k[0])+"/CONTCAR")
                        #         i.bundle.append(i.detail["path"]+"/cc_"+str(j[0])+"_"+str(k[0]))
                        i.state = "done"

                elif i.detail["type"] == "distort":
                    if os.path.isdir(i.detail["path"]):
                        i.state = "done"
                    elif i.detail["structure"] in self.job_path_list and self.find_job_by_path(
                            i.detail["structure"]).state == "done":
                        p = poscar.poscar()
                        if self.find_job_by_path(i.detail["structure"]).bundle:
                            pass
                        else:
                            p.read(i.detail["structure"] + "/CONTCAR")
                            scc_lib.linux_command("mkdir " + i.detail["path"])
                            p.distort(i.detail["mode"]).write(i.detail["path"] + "/CONTCAR")
                        i.state = "done"

                elif i.detail["type"] == "distort_sphere":
                    if os.path.isdir(i.detail["path"]):
                        i.state = "done"
                    elif (i.detail["structure"] in self.job_path_list
                          and self.find_job_by_path(i.detail["structure"]).state == "done"
                          and (not "/" in i.detail["distort_list"]
                               or ("/" in i.detail["distort_list"]
                                   and self.find_job_by_path(i.detail["distort_list"]).state == "done"
                                   and self.find_job_by_path(origin_name(i.detail["distort_list"])).state == "done"))):

                        if i.detail["distort_list"] == "0":
                            i.state = "fail"
                        else:
                            if "/" in i.detail["distort_list"]:
                                distort_poscar = poscar.poscar()
                                origin_poscar = poscar.poscar()
                                origin_poscar.read(origin_name(i.detail["distort_list"]) + "/CONTCAR")
                                distort_poscar.read(i.detail["distort_list"] + "/CONTCAR")
                                if "half" in i.detail["path"]:
                                    distortion = (distort_poscar - origin_poscar) * 0.5
                                elif "double" in i.detail["path"]:
                                    distortion = (distort_poscar - origin_poscar) * 2
                                else:
                                    distortion = distort_poscar - origin_poscar
                                distort_list = []
                                current_element = False
                                neis = 0
                                for j in distort_poscar.distance(distort_poscar.label[0]):
                                    if not current_element == j[0][0]:
                                        neis += 1
                                        if neis > 3:
                                            break
                                        else:
                                            current_element = j[0][0]
                                    abs_origin = (origin_poscar.position[distort_poscar.label.index(j[0])]
                                                  - origin_poscar.position[0]).dot(origin_poscar.axis)
                                    abs_distortion = distortion.position[distort_poscar.label.index(j[0])].dot(
                                        distortion.axis)
                                    if np.linalg.norm(abs_origin):
                                        distort_list.append(
                                            round(abs_origin.dot(abs_distortion) / np.linalg.norm(abs_origin), 3))
                            else:
                                distort_list = [float(j) for j in i.detail["distort_list"].split("_")]
                            high_sym_poscar = poscar.poscar()
                            high_sym_poscar.read(i.detail["structure"] + "/CONTCAR")
                            high_sym_distance = high_sym_poscar.distance(high_sym_poscar.label[0])
                            if (distort_list[5] - distort_list[0]) < 1e-5:
                                distort_mode = "oh"
                            elif ((distort_list[3] - distort_list[0]) < 1e-3
                                  and (distort_list[5] - distort_list[4]) < 1e-3):
                                distort_mode = "ez"
                            elif ((distort_list[1] - distort_list[0]) < 1e-3
                                  and (distort_list[5] - distort_list[2]) < 1e-3):
                                distort_mode = "cz"
                            else:
                                distort_mode = "nosym"

                            avg = 0
                            z_list = []
                            for j in range(7)[1:]:
                                avg += high_sym_distance[j][1]
                                pos = (high_sym_poscar.position[high_sym_poscar.label.index(high_sym_distance[j][0])]
                                       - high_sym_poscar.position[0])
                                if abs(pos[0]) < 1e-3 and abs(pos[1]) < 1e-3:
                                    z_list.append(j)
                            avg /= 6
                            pre_distort_list = []
                            if distort_mode == "ez":
                                for j in range(7)[1:]:
                                    pre_distort_list.append(
                                        avg - high_sym_distance[j][1] + (1e-5 if j in z_list else -1e-5))
                            elif distort_mode == "cz":
                                for j in range(7)[1:]:
                                    pre_distort_list.append(
                                        avg - high_sym_distance[j][1] - (1e-5 if j in z_list else -1e-5))
                            scc_lib.linux_command("mkdir " + i.detail["path"])
                            high_sym_poscar.distort_sphere(pre_distort_list).distort_sphere(distort_list).write(
                                i.detail["path"] + "/CONTCAR")
                            i.state = "done"

                elif i.detail["type"] == "expand":
                    if os.path.isdir(i.detail["path"]):
                        i.state = "done"
                    elif self.find_job_by_path(i.detail["structure"]).state == "done":
                        p = poscar.poscar()
                        if self.find_job_by_path(i.detail["structure"]).bundle:
                            pass
                        else:
                            p.read(i.detail["structure"] + "/CONTCAR")
                            scc_lib.linux_command("mkdir " + i.detail["path"])
                            p.expand(i.detail["tm"])
                            p.write(i.detail["path"] + "/CONTCAR")
                        i.state = "done"

                elif "slab" in i.detail["type"]:
                    if os.path.isdir(i.detail["path"]):
                        i.state = "done"
                    elif self.find_job_by_path(i.detail["structure"]).state == "done":
                        p = poscar.poscar()
                        if self.find_job_by_path(i.detail["structure"]).bundle:
                            pass
                        else:
                            p.read(i.detail["structure"] + "/CONTCAR")
                            scc_lib.linux_command("mkdir " + i.detail["path"])
                            if i.detail["type"] == "slab":
                                p.slab(i.detail["vac_length"])
                            elif i.detail["type"] == "slab2":
                                p.slab2(i.detail["vac_length"])
                            p.write(i.detail["path"] + "/CONTCAR")
                        i.state = "done"

                elif "primcell" in i.detail["type"]:
                    if os.path.isdir(i.detail["path"]):
                        i.state = "done"
                    elif self.find_job_by_path(i.detail["structure"]).state == "done":
                        scc_lib.linux_command("mkdir " + i.detail["path"])
                        scc_lib.linux_command(
                            "cp " + i.detail["structure"] + "/CONTCAR " + i.detail["path"] + "/POSCAR")
                        cwd = os.getcwd()
                        os.chdir(i.detail["path"])
                        scc_lib.linux_command("vaspkit -task 303")
                        os.chdir(cwd)
                        scc_lib.linux_command(
                            "mv " + i.detail["path"] + "/PRIMCELL.vasp " + i.detail["path"] + "/CONTCAR")
                        i.state = "done"

                elif i.detail["type"] == "plot_pot":
                    if os.path.isfile(i.detail["path"] + ".png"):
                        i.state = "done"
                    elif self.find_job_by_path(i.detail["pot"]).state == "done":
                        cwd = os.getcwd()
                        os.chdir(i.detail["pot"])
                        scc_lib.linux_command("(echo 426;echo 3)|vaspkit")
                        os.chdir(cwd)

                        if "atom" in i.detail["path"]:
                            ppot.atom_pot(i.detail["pot"], i.detail['path'])
                        else:
                            procar = read_procar(i.detail["pot"])[0]
                            for j in procar:
                                if j["occ"] < 0.3:
                                    vbm = procar[procar.index(j) - 1]["energy"]
                                    break
                            # file.write("\n"+i.detail["path"]+" "+str(vbm)+" ")

                            with open(i.detail["pot"] + "/PLANAR_AVERAGE.dat", "r") as f:
                                lines = f.readlines()
                            match_result = re.search("#(.+), (.+)", lines[0])
                            if not match_result:
                                match_result = re.search("#(.+) (.+)", lines[0])
                            xlabel, ylabel = match_result.groups()
                            x = [xlabel]
                            y = [ylabel]
                            for j in lines[1:]:
                                count = 0
                                for k in j.replace("\n", "").split(" "):
                                    if k:
                                        count += 1
                                        if count % 2 == 1:
                                            x.append(float(k))
                                        else:
                                            y.append(float(k))

                            count = {}
                            for j in [round(j, 3) for j in y[1:]]:
                                if j in count:
                                    count[j] += 1
                                else:
                                    count[j] = 1

                            count = 1
                            ysum = y[1]
                            temp_count = []
                            temp_ysum = []
                            file.write("avg " + i.detail["pot"] + "\n")
                            for j in range(len(x))[2:-1]:
                                count += 1
                                ysum += y[j]
                                if (y[j] - y[j - 1]) / (x[j] - x[j - 1]) * (y[j + 1] - y[j]) / (x[j + 1] - x[j]) < 0:
                                    temp_count.append(count)
                                    temp_ysum.append(ysum)
                                    # file.write(str(count)+" "+str(ysum)+" ")
                                    count = 0
                                    ysum = 0

                            for j in temp_count:
                                file.write(str(j) + " ")
                            file.write("\n")
                            for j in temp_ysum:
                                file.write(str(j) + " ")
                            file.write("\n")
                            for j, k in zip(temp_count, temp_ysum):
                                file.write(str(k / j) + " ")
                            file.write("\n")

                            count = 0
                            ysum = 0
                            file.write("sum " + i.detail["pot"] + " ")
                            for j in range(len(x))[1:]:
                                count += 1
                                ysum += y[j]
                                file.write(str(ysum / count) + " ")
                            file.write("\n")

                            good_range = 3
                            good_list = []
                            for j in range(len(x))[good_range + 1:len(x) - good_range]:
                                if abs((y[j] - y[j - good_range]) / (x[j] - x[j - good_range]) -
                                       (y[j] - y[j + good_range]) / (x[j] - x[j + good_range])) < 0.01:
                                    good_list.append([x[j], y[j]])
                            if len(good_list) > 5:
                                file.write("pot " +
                                           i.detail["pot"] + " " +
                                           str(sum(list(zip(*good_list))[1]) / len(good_list)) + " " +
                                           str(good_list[int(len(good_list) / 2)][1]) + "\n")
                            else:
                                # file.write("pot "+
                                #            i.detail["pot"]+" "+
                                #            str(y[1])+"\n")
                                pass

                            plot_pot(x, y, i.detail["path"])
                        i.state = "done"

                elif i.detail["type"] == "plot_cc":
                    if not i.bundle:
                        # vars in i.detail
                        for j in ["grd", "ex", "ex_at_grd", "grd_at_ex"]:
                            i.bundle.append(self.find_job_by_path(i.detail[j]))
                    if set([j.state for j in i.bundle]) == {"done"}:
                        file.write("\n" + i.detail["path"] + "\n")

                        energy_list = [energy(j.detail["path"]) if j else j for j in i.bundle]
                        if not os.path.isfile(i.detail["path"] + ".png"):
                            plot_cc(energy_list, i.detail["path"])
                        i.state = "done"
                elif i.detail["type"] == "plot_ctl":
                    if not i.bundle:
                        ctl_pattern = re.compile(i.detail["pattern"].replace("xxx", "(\w)(\d)") + "$")
                        for j in self.job_list:
                            if "path" in j.detail:
                                match_result = ctl_pattern.search(j.detail["path"])
                                if match_result:
                                    i.bundle.append(j)
                        i.bundle.append(self.find_job_by_path(i.detail["pattern"].replace("xxx", "")))
                        i.bundle.append(self.find_job_by_path(i.detail["unitcell"]))
                    # if set([j.state for j in i.bundle])=={"done"}:
                    if i.bundle[-1].state == "done":
                        # bandgap
                        job_unitcell = self.find_job_by_path(i.detail["unitcell"])
                        procar = read_procar(job_unitcell.detail["path"])[0]
                        for j in procar:
                            if j["occ"] < 0.3:
                                cbm = j["energy"]
                                vbm = procar[procar.index(j) - 1]["energy"]
                                break
                        energy_unitcell = energy(job_unitcell.detail["path"])
                        # 是超胞
                        if "supercell" in self.find_job_by_path(i.detail["path"].split("/")[0] + "/sub").detail[
                            "structure"]:
                            expand_job = self.find_job_by_path(
                                "0" + i.detail["path"].split("/")[0].split("_")[0][1:] + "/supercell")
                            energy_unitcell *= np.linalg.det(expand_job.detail["tm"])
                        energy_list = [[-j.detail["nele"], energy(j.detail["path"]) - energy_unitcell] for j in
                                       i.bundle[:-1] if j.state == "done"]

                        energy_list.sort(key=lambda x: x[0])
                        file.write("\n[" + str(vbm) + "," + str(cbm) + ",")
                        for j in range(len(energy_list) - 1):
                            file.write(str(energy_list[j][1] - energy_list[j + 1][1]) + ",")
                        file.write("\"" + i.detail["path"] + "\"],\n")
                        for j in energy_list:
                            file.write(str(j[0]) + " " + str(j[1]) + "\n")
                        if not os.path.isfile(i.detail["path"] + ".png"):
                            plot_ctl([vbm, cbm], energy_list, i.detail["path"])
                        i.state = "done"

                elif i.detail["type"] == "plot_ctl3":
                    if not i.bundle:
                        for j in [i.detail["p1"],
                                  i.detail["p1"].replace("p1", ""),
                                  i.detail["m1"],
                                  i.detail["m1"].replace("m1", ""),
                                  i.detail["unitcell"]]:
                            i.bundle.append(self.find_job_by_path(j))

                    if set([j.state for j in i.bundle]) in [{"done", "wrong"}, {"done"}]:
                        # bandgap
                        procar = read_procar(i.bundle[-1].detail["path"])[0]
                        for j in procar:
                            if j["occ"] < 0.3:
                                cbm = j["energy"]
                                vbm = procar[procar.index(j) - 1]["energy"]
                                break

                        energy_list = [[-i.bundle[0].detail["nele"],
                                        energy(i.bundle[0].detail["path"]) - energy(i.bundle[1].detail["path"])],
                                       [-i.bundle[1].detail["nele"], 0],
                                       [-i.bundle[2].detail["nele"],
                                        energy(i.bundle[2].detail["path"]) - energy(i.bundle[3].detail["path"])]]

                        energy_list.sort(key=lambda x: x[0])
                        file.write("\n[" + str(vbm) + "," + str(cbm) + ",")
                        for j in range(len(energy_list) - 1):
                            if i.bundle[j].state == "wrong":
                                file.write(str(1000 + energy_list[j][1] - energy_list[j + 1][1]) + ",")
                            else:
                                file.write(str(energy_list[j][1] - energy_list[j + 1][1]) + ",")
                        file.write("\"" + i.detail["path"] + "\"],\n")
                        for j in energy_list:
                            file.write(str(j[0]) + " " + str(j[1]) + "\n")
                        if not os.path.isfile(i.detail["path"] + ".png"):
                            plot_ctl([vbm, cbm], energy_list, i.detail["path"])
                        i.state = "done"

                elif i.detail["type"] == "ctl4" and not delay:

                    if self.find_job_by_path(
                            i.detail["system"].replace("1", "0").split("_")[0] + "/unitcell_" + i.detail[
                                "level"]).state == "done":
                        energy_dict1 = {}
                        energy_dict2 = {}
                        for j in self.job_path_list:
                            for k in ["grdm1", "grdp1", "ex", "grd"]:
                                if i.detail["system"] + "/relax_" + k in j and i.detail["level"] in j:
                                    if self.find_job_by_path(j).state == "done":
                                        energy_dict1[j.replace(i.detail["system"] + "/relax_", "").replace(
                                            "_" + i.detail["level"], "")] = energy(j) - charge_corr_dict[
                                            i.detail["system"].replace("1", "0").split("_")[0]] * \
                                                                            self.find_job_by_path(j).detail["nele"] ** 2
                                    elif self.find_job_by_path(j).state == "wrong":
                                        energy_dict1[j.replace(i.detail["system"] + "/relax_", "").replace(
                                            "_" + i.detail["level"], "")] = energy(j) + 1e5 - charge_corr_dict[
                                            i.detail["system"].replace("1", "0").split("_")[0]] * \
                                                                            self.find_job_by_path(j).detail["nele"] ** 2
                                    break

                        if "grd" in energy_dict1:
                            for j in energy_dict1:
                                if "grdp1" in j:
                                    energy_dict2["0/1" + j.replace("grdp1", "")] = energy_dict1[j] - energy_dict1["grd"]

                        if "grdm1" in energy_dict1:
                            for j in energy_dict1:
                                if not "grdp1" in j and not "grdm1" in j:
                                    energy_dict2["-1/0" + j.replace("grd", "").replace("ex", "")] = energy_dict1[j] - \
                                                                                                    energy_dict1[
                                                                                                        "grdm1"]

                        procar = read_procar(
                            i.detail["system"].replace("1", "0").split("_")[0] + "/unitcell_" + i.detail["level"])[0]
                        for j in procar:
                            if j["occ"] < 0.3:
                                energy_dict2["CBM"] = j["energy"]
                                energy_dict2["VBM"] = procar[procar.index(j) - 1]["energy"]
                                break

                        energy_dict2["system"] = i.detail["system"]

                        file.write(str(energy_dict2) + "ctl4\n")
                        i.state = "done"

                        # update plot_nei要在所有bundle的vasp计算之后
                elif i.detail["type"] == "plot_nei" and not delay:
                    for j in self.job_list:
                        if ("path" in j.detail
                                and j.detail["path"].split("/")[0] == i.detail["path"].split("/")[0]
                                and "relax" in j.detail["path"]
                                and j.state == "done"
                                and not j.bundle):
                            i.bundle.append(j)
                    # 多于一个才可以画子图
                    if len(i.bundle) > 0:
                        file.write("\n" + i.detail["path"] + "\n")
                        p = poscar.poscar()
                        nei_list = []
                        match_result = re.search("(\w+)_(\d+)", i.detail["atom"])
                        for j in i.bundle:
                            p.read(j.detail["path"] + "/CONTCAR")
                            current_nei = []
                            for k in p.env((match_result.group(1), int(match_result.group(2))), 2, False):
                                if k[1] < i.detail["range"]:
                                    current_nei.append(k)
                                else:
                                    current_nei.append(j.detail["path"])
                                    break
                            nei_list.append(current_nei)
                        for j in nei_list:
                            file.write(str(j) + ",\n")
                        if len(i.bundle) > 1:
                            plot_nei(nei_list, i.detail["path"])
                        i.state = "done"

                elif i.detail["type"] == "plot_band":
                    if os.path.isfile(i.detail["path"] + ".png"):
                        i.state = "done"
                    elif self.find_job_by_path(i.detail["structure"]).state == "done":
                        cwd = os.getcwd()
                        os.chdir(i.detail["structure"])
                        scc_lib.linux_command("vaspkit -task 211")
                        os.chdir(cwd)
                        band.plot_band(i.detail["structure"] + "/BAND.dat", i.detail["structure"] + "/KLABELS",
                                       i.detail["path"] + ".png")
                        i.state == "done"

        return None

    def run(self, file):

        for i in list(set([i.split("/")[0] for i in self.job_path_list])):
            if not os.path.isdir(i):
                scc_lib.linux_command("mkdir " + i)

        # 清理过时任务

        # jobid_template=re.compile("(\d+)")
        # for i in os.listdir():
        #     if os.path.isdir(i):
        #         for j in os.listdir(i):
        #             if os.path.isdir(i+"/"+j):
        #                 listdir=os.listdir(i+"/"+j)
        #                 if "CONTCAR" in listdir and "PROCAR" in listdir and "jobid" in listdir and i[0] in "123456789" and "OSZICAR" in listdir:
        #                     with open(i+"/"+j+"/jobid","r") as f:
        #                         jobid=jobid_template.search(f.readlines()[-1]).group(1)
        #                         print(i+"/"+j,end=" ")
        #                     if jobid+".log" in listdir:
        #                         if "relax" in j:
        #                             with open(i+"/"+j+"/"+jobid+".log","r") as f:
        #                                 log_lines=f.readlines()
        #                             if "for stderr output of this job" in log_lines[-2]:
        #                                 log_lines=log_lines[:-6]
        #                             if "reached required accuracy - stopping structural energy minimisation" in log_lines[-2]:
        #                                 p=poscar.poscar()
        #                                 p.read(i+"/"+j+"/CONTCAR")
        #                                 print(p.distance(p.label[0],False)[:7],end=" ")
        #                                 print(energy(i+"/"+j),end=" ")
        #                                 procar=read_procar(i+"/"+j)
        #                                 spd1=(i.split("_")[1],"s")
        #                                 for k in procar[-1]:
        #                                     if k["occ"]<0.7:
        #                                         if not spd1 in k["spd"]:
        #                                             print("ws",end=" ")
        #                                         elif k["occ"]>0.3:
        #                                             print("ws",end=" ")
        #                                         break

        #                                 spd2=spd1=(i.split("_")[1],"p")
        #                                 for k in procar[0][::-1]:
        #                                     if k["occ"]>0.3:
        #                                         if not spd2 in k["spd"]:
        #                                             print("wp",end=" ")
        #                                         elif k["occ"]<0.7:
        #                                             print("wp",end=" ")
        #                                         break
        #                                 print()
        #                             else:
        #                                 print("notreach")
        #                         else:
        #                             p=poscar.poscar()
        #                             p.read(i+"/"+j+"/CONTCAR")
        #                             print(p.distance(p.label[0],False)[:7],end=" ")
        #                             print(energy(i+"/"+j),end=" ")
        #                             procar=read_procar(i+"/"+j)
        #                             spd1=(i.split("_")[1],"s")
        #                             for k in procar[-1]:
        #                                 if k["occ"]<0.7:
        #                                     if not spd1 in k["spd"]:
        #                                         print("ws",end=" ")
        #                                     elif k["occ"]>0.3:
        #                                         print("ws",end=" ")
        #                                     break

        #                             spd2=spd1=(i.split("_")[1],"p")
        #                             for k in procar[0][::-1]:
        #                                 if k["occ"]>0.3:
        #                                     if not spd2 in k["spd"]:
        #                                         print("wp",end=" ")
        #                                     elif k["occ"]<0.7:
        #                                         print("wp",end=" ")
        #                                     break
        #                             print()
        #                     else:
        #                         print("nolog")

        # if not i+"/"+j in self.job_path_list:
        #     if not "notinlist" in i+"/"+j:
        #         if "jobid" in listdir:
        #             with open(i+"/"+j+"/jobid","r") as f:
        #                 scc_lib.linux_command("bkill "+jobid_template.search(f.readlines()[-1]).group(1))
        #         scc_lib.linux_command("mv "+i+"/"+j+" "+i+"/"+j+"_nil")
        #     file.write(i+"/"+j+"\n")

        self.update(file, True)
        gonna_finish = False
        while True:
            job_left = {}
            job_left["wait"] = []
            job_left["run"] = []
            job_left["done"] = []
            job_left["fail"] = []
            job_left["wrong"] = []
            job_dict_dict = {}

            self.update(file, False)
            for i in self.job_list:
                if "path" in i.detail:
                    job_left[i.state].append(i.detail["path"])
                    if i.detail["path"].split("/")[0][0] in ["1", "2"] and i.detail["type"] == "vasp":
                        if not i.detail["path"].split("/")[0] in job_dict_dict:
                            job_dict_dict[i.detail["path"].split("/")[0]] = {}
                            job_dict_dict[i.detail["path"].split("/")[0]]["system"] = i.detail["path"].split("/")[0]
                        job_dict_dict[i.detail["path"].split("/")[0]][i.detail["path"].split("/")[1]] = i.state

            # with open("state", "w") as f:
            #     print(job_dict_dict.values())
            #     f.write(str(table(list(job_dict_dict.values()))) + "\n")

            # if not job_left["run"]:
            #     if gonna_finish:
            #         break
            #     else:
            #         gonna_finish = True
            #         for i in ["run", "fail", "wrong"]:
            #             print(str(len(job_left[i])) + " " + i)
            #             for j in job_left[i]:
            #                 print(j)
            #             print()
            # else:
            #     gonna_finish = False
            #     for i in ["run", "fail", "wrong"]:
            #         print(str(len(job_left[i])) + " " + i)
            #         for j in job_left[i]:
            #             print(j)
            #         print()
            time.sleep(600)
        return None


def table(dict_list):
    key = []
    table = []
    for i in dict_list:
        for j in i.keys():
            if not j in key:
                key.append(j)
    for i in dict_list:
        for j in key:
            if not j in i:
                i[j] = None
    table.append(key)
    for i in dict_list:
        table.append([i[j] for j in key])

    for i in table:
        for j in i:
            print(j, end=" ")
        print()
    return None


if __name__ == '__main__':
    with open("report", "w") as f:
        jobs().run(f)
