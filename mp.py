#!/home/phys/qif/anaconda3/bin/python

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:38:51 2021

@author: dugue
"""

import sys
from pymatgen.ext.matproj import MPRester

# ['band_gap', 'cif', 'density', 'diel', 'e_above_hull', 'elasticity',
#  'elements', 'energy', 'energy_per_atom',
#  'formation_energy_per_atom', 'full_formula', 'hubbards',
#  'icsd_id', 'icsd_ids', 'is_compatible', 'is_hubbard',
#  'material_id', 'nelements', 'nsites', 'oxide_type', 'piezo',
#  'pretty_formula', 'spacegroup', 'tags', 'task_ids',
#  'total_magnetization', 'unit_cell_formula', 'volume']
prop_list = ['material_id',
             'pretty_formula',
             'energy_per_atom',
             'formation_energy_per_atom',
             'band_gap',
             'symbol',
             'point_group',
             'total_magnetization']

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


# 交互/非交互访问material project
def matproj(formula):
    with MPRester("VUswQqBEWe4VFBZD25") as m:
        prop_dict = []
        for i in m.get_data(formula):
            prop_dict.append({})
            for j in i:
                if j == "spacegroup":
                    prop_dict[-1]["point_group"] = point_group_dict[i["spacegroup"]["point_group"]]
                    prop_dict[-1]["symbol"] = i["spacegroup"]["symbol"]
                elif j in prop_list:
                    if j == "material_id":
                        m.get_structure_by_material_id(i[j],
                                                       final=True,
                                                       conventional_unit_cell=True).to(
                            filename=formula + i[j].replace("mp", ""),
                            fmt="poscar")
                    prop_dict[-1][j] = i[j]

    prop_dict.sort(key=lambda x: x['formation_energy_per_atom'])

    for i in prop_list:
        print("%-10s" % i[:10], end=" ")
    print()

    for i in prop_dict:
        for j in prop_list:
            if type(i[j]) == float:
                i[j] = round(i[j], 5)
            print("%-10s" % i[j], end=" ")
        print()


if __name__ == '__main__':
    # matproj("Eu-O")
    matproj("O")

    # from mp_api.client import MPRester
    # with MPRester("mCQcJ91fLvg014UNkWYP5BDQOEp3CK3U") as mpr:
    #     docs = mpr.summary.search(material_ids=["mp-2605"])

    # from mp_api.client import MPRester

    # with MPRester("mCQcJ91fLvg014UNkWYP5BDQOEp3CK3U") as mpr:
    #     docs = mpr.summary.search(material_ids=["mp-149"], fields=["structure"])
    #     structure = docs[0].structure
    #     # -- Shortcut for a single Materials Project ID:
    #     structure = mpr.get_structure_by_material_id("mp-149")

    # from mp_api.client import MPRester
    #
    # with MPRester("mCQcJ91fLvg014UNkWYP5BDQOEp3CK3U") as mpr:
    #     ph_bs = mpr.get_phonon_bandstructure_by_material_id("mp-149")

    # curl - X
    # 'GET' \
    # 'https://api.materialsproject.org/materials/summary/mp-2605/?_all_fields=false' \
    # - H
    # 'accept: application/json' \
    # - H
    # 'X-API-KEY: mCQcJ91fLvg014UNkWYP5BDQOEp3CK3U'
    #
    # import requests
    #
    # url = "https://api.materialsproject.org/materials/summary/mp-2605/?_all_fields=false"
    # headers = {
    #     "accept": "application/json",
    #     "X-API-KEY": "mCQcJ91fLvg014UNkWYP5BDQOEp3CK3U"
    # }
    #
    # response = requests.get(url, headers=headers)
    # print(response.json())
