import json
import pymatgen.core as mg

with open("mpb","r") as f:
    while True:
        i=f.readline()
        if "search" in i:
            continue
        else:
            j=json.loads(i)
            print(j["material_id"],j["material_id"],j["pretty_formula"],j["formation_energy_per_atom"])
            structure = mg.Structure.from_str(j["cif"], fmt="cif")
            structure.to(filename=f"{j['material_id']}.vasp",fmt="poscar")