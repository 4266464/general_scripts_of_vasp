import json
import os
import re

import numpy as np
from pymatgen.ext.matproj import MPRester

point_group_dict = {"1": "C1", "4": "C4", "2": "C2", "2/m": "C2h", "mm2": "C2v", "3": "C3", "-6": "C3h", "-3": "C3i",
                    "3m": "C3v", "4/m": "C4h", "4mm": "C4v", "6": "C6", "6/m": "C6h", "6mm": "C6v", "-1": "Ci",
                    "m": "Cs", "222": "D2", "-42m": "D2d", "mmm": "D2h", "32": "D3", "-3m": "D3d", "-6m2": "D3h",
                    "422": "D4", "4/mmm": "D4h", "622": "D6", "6/mmm": "D6h", "432": "O", "m-3m": "Oh", "-4": "S4",
                    "23": "T", "-43m": "Td", "m-3": "Th"}

properties_displayed = ['formation_energy_per_atom', 'pretty_formula', 'material_id']


def convert_formula_to_table(formula):
    number_of_element = len(list(filter(str.isupper, formula)))
    expression = "^" + number_of_element * "([A-Z][a-z]?)(\d?)" + "$"
    print(expression, formula)
    match_group = re.match(expression, formula).groups()
    match_group = np.array(match_group).reshape(-1, 2)
    match_group[..., 1] = np.where(match_group[..., 1] == "", 1, match_group[..., 1]).astype(int)

    dtype = np.dtype([("element", "U2"), ("number", "<i4")])
    table = np.array([("El", 0)] * number_of_element, dtype=dtype)
    table["element"] = match_group[..., 0]
    table["number"] = match_group[..., 1]
    table = np.sort(table, order="element")

    # print(f"formula {formula} has been convert to {table}.")
    return table


class material_project_repo:
    def __init__(self, repo_path):
        self.path = repo_path
        self.set_list = [name.split(".")[0] for name in os.listdir(self.path)]

    def check(self):
        '''
        Check whether the material project data can be loaded by JSON.
        '''
        for set in self.set_list:
            try:
                json.load(open(self.path + "/" + set + ".json", "r"))
            except json.decoder.JSONDecodeError:
                print(set)

    def query(self, compound):
        assert isinstance(compound, str), f"compound {compound} should be a string."

        # the compound in form of element set or formula
        is_set = "-" in compound
        if is_set:
            current_set = "-".join(sorted(compound.split("-")))
        else:
            table = convert_formula_to_table(compound)
            current_set = "-".join(table["element"])

        if current_set not in self.set_list:
            print(f"compound {compound} not found in local repo {self.path}. updating...")
            if not self.update(current_set):
                return None

        result_list = json.load(open(self.path + "/" + current_set + ".json", "r"))

        def print_result(result_list):
            print(properties_displayed)
            for result in result_list:
                print([result[tag] for tag in result if tag in properties_displayed])

        if is_set:
            print_result(result_list)
            return result_list
        else:
            def result_filter(result):
                return np.all(convert_formula_to_table(result["pretty_formula"])["number"] == table["number"])

            result_list = list(filter(result_filter, result_list))
            if result_list:
                print_result(result_list)
                return result_list
            else:
                print(f"homologous series are found instead of the compound {compound}, "
                      "which indicates the compound is not included in the material project database. You could try quering online.")
                return None

    def update(self, compound_set):
        '''
        update all the compound in the repo
        :param self.path: where the repo is
        :return: None
        '''
        with MPRester("VUswQqBEWe4VFBZD25") as m:
            result_list = m.get_data(compound_set)
            if result_list:
                json.dump(result_list, open(self.path + "/" + compound_set + ".json", "w"))
                return True
            else:
                return False


if __name__ == "__main__":
    # table = convert_formula_to_table("Cs2NaBiCl6")
    # print(table, table["element"], table["number"])
    mp = material_project_repo("material_project_repo")
    # mp.query("AgBr")
    for i in "CaWO4 SrWO4 BaWO4".split(" "):
        print(i)
        mp.query(i)

# TODO: 没有的也存下来， 不然每次都找