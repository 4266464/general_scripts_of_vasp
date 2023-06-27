from old_poscar import poscar
import re
from match_times import match_times

def read_procar(path):
    band_pattern = re.compile("band\s+(\d+) # energy\s+([\d\.\-]+) # occ.\s+([\d\.]+)")
    p = poscar()
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