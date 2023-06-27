# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 23:46:31 2021

@author: dugue
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_ctl(charged_energy_dict, temperature, path, band1=[-10, 20, 100], fermi_where=False):
    # 使用价态而不是电子数
    y = {i: [] for i in charged_energy_dict}
    fermi_level = False
    x = np.linspace(*band1)
    for i in x:
        total_charge = 0
        if not fermi_level:
            for j in charged_energy_dict:
                y[j].append(min([k[0] * i + k[1] for k in charged_energy_dict[j]]))
                total_charge += sum(
                    [k[0] * np.exp(-(k[0] * i + k[1]) / (8.617e-05 * temperature)) for k in charged_energy_dict[j]])
            if total_charge < 0:
                fermi_level = i
        else:
            for j in charged_energy_dict:
                y[j].append(min([k[0] * i + k[1] for k in charged_energy_dict[j]]))

    y_max = max(max(i) for i in y.values())
    y_min = min(min(i) for i in y.values())

    latex = {"vcs": "v\-(Cs)",
             "vna": "v\-(Na)",
             "vbi": "v\-(Bi)",
             "vcl": "v\-(Cl)",
             "bics": "Bi\-(Cs)",
             "nacs": "Na\-(Cs)",
             "csna": "Cs\-(Na)",
             "bina": "Bi\-(Na)",
             "csbi": "Cs\-(Bi)",
             "nabi": "Na\-(Bi)",
             "icl": "i\-(Cl)",
             "ina": "i\-(Na)",
             "ics": "i\-(Cs)",
             "AlBa": "Al\-(Ba)",
             "AlCa": "Al\-(Ca)",
             "AlSr": "Al\-(Sr)",
             "BaAl": "Ba\-(Al)",
             "BiAl": "Bi\-(Al)",
             "BiBa": "Bi\-(Ba)",
             "BiCa": "Bi\-(Ca)",
             "BiSr": "Bi\-(Sr)",
             "CaAl": "Ca\-(Al)",
             "FeAl": "Fe\-(Al)",
             "FeCa": "Fe\-(Ca)",
             "MnSr": "Mn\-(Sr)",
             "PbAl": "Pb\-(Al)",
             "PbCa": "Pb\-(Ca)",
             "SrAl": "Sr\-(Al)",
             "vAl": "v\-(Al)",
             "vBa": "v\-(Ba)",
             "vCa": "v\-(Ca)",
             "vO": "v\-(O)",
             "vSr": "v\-(Sr)",
             "iAl": "i\-(Al)",
             "iBi": "Bi\-(Ca)+i\-(O)",
             "iCa": "i\-(Ca)",
             "iO": "i\-(O)",
             "Bicca": "Bi\-(Ca)+Ca\-(Al)",
             "Bicva": "Bi\-(Ca)+V\-(Al)",
             "Fe5Al": "Fe\-(Al\-(5))",
             "Fe4Al": "Fe\-(Al\-(4))",
             "Fe6Al": "Fe\-(Al\-(6))",
             "Mn6Al": "Mn\-(Al\-(6))",
             "Mn4Al": "Mn\-(Al\-(4))",
             "BiZr": "Bi\-(Zr)",
             "BivZrCl": "Bi\-(Zr)+v\-(Cl)",
             "i2Cl": "i\-(Cl)(Cl\-(2))",
             "iCl": "i\-(Cl)(Cl\-(6))",
             "iCs": "i\-(Cs)",
             "SbZr": "Sb\-(Zr)",
             "SbvZrCl": "Sb\-(Zr)+v\-(Cl)",
             "SeZr": "Se\-(Zr)",
             "TeZr": "Te\-(Zr)",
             "vCl": "v\-(Cl)",
             "vCs": "v\-(Cs)",
             }

    lines = ["", "", "", ""]
    mat = []
    for i in y:
        temp_mat = []
        # plt.plot(x, y[i], label=latex[i])
        # lines[0]+=f"Fermi_energy {latex[i]}  "
        lines[0] += f"Fermi_energy {i}  "
        lines[1] += "eV eV  "
        temp_mat.append([0, f"{0} {y[i][0]:.2f}  "])
        temp_mat.append([x[-1] - x[0], f"{x[-1] - x[0]:.2f} {y[i][-1]:.2f}  "])
        mat.append(temp_mat)
        # lines[2]+=f"{0} {y[i][0]:.2f}  "
        # lines[3] += f"{x[-1]-x[0]:.2f} {y[i][-1]:.2f}  "
        # print(x[0],y[i][0],x[-1],y[i][-1],i)
        plt.plot(x, y[i], label=i)
        # plt.annotate(i, xy=(x[0], y[i][0]), xytext=(0, 0), textcoords='offset points', fontsize=6)
    vbm = x[0]

    crossing = {}
    for i in charged_energy_dict:
        crossing[i] = {}
        charged_energy_dict[i].sort(key=lambda x: x[0])
        for j in range(len(charged_energy_dict[i]) - 1):
            print(charged_energy_dict[i])
            temp_x = -(charged_energy_dict[i][j][1] - charged_energy_dict[i][j + 1][1]) / \
                     (charged_energy_dict[i][j][0] - charged_energy_dict[i][j + 1][0])
            temp_y = temp_x * charged_energy_dict[i][j][0] + charged_energy_dict[i][j][1]
            crossing[i][f"{charged_energy_dict[i][j][0]}/{charged_energy_dict[i][j + 1][0]}"] = [temp_x, temp_y]

    for i in range(max([len(crossing[i]) for i in crossing])):
        lines.append("")
    # print(crossing,mat)
    for i, j in zip(crossing, range(len(mat))):
        for k in crossing[i]:
            mat[j].append([crossing[i][k][0] - vbm,
                           f"{crossing[i][k][0] - vbm:.2f} {crossing[i][k][1]:.2f} {k.replace('/', '|')} "])
            # print(f"{i} {j} {crossing[i][j][0]} {crossing[i][j][1]}")
            # plt.plot([crossing[i][j], crossing[i][j]], [y_min * (1 + factor), y_max * (1 - factor)], color="blue",
            # linestyle="dashed")
            # plt.scatter(crossing[i][j][0],crossing[i][j][1])
            plt.annotate(f"{k}", xy=(crossing[i][k][0], crossing[i][k][1]), xytext=(-5, 5), textcoords='offset points',
                         fontsize=8)
        for k in range(len(lines) - 4 - len(crossing[i])):
            # lines[j]+=3*" "
            mat[j].append([10000, 3 * " "])

    for i in range(len(mat)):
        mat[i].sort(key=lambda x: x[0])

    for i in mat:
        for j, k in zip(i, range(len(lines))[2:]):
            lines[k] += j[1]

    print(20 * "-")
    for i in lines:
        print(i)
    print(20 * "-")

    factor = 2e-4
    # if band2:
    #     plt.plot([band2[0], band2[0]], [y_min * (1 + factor), y_max * (1 - factor)], color="blue")
    #     plt.plot([band2[1], band2[1]], [y_min * (1 + factor), y_max * (1 - factor)], color="blue")
    #     plt.annotate("VBM2", xy=(band2[0], y_min), xytext=(0, 0), textcoords='offset points', fontsize=6)
    #     plt.annotate("CBM2", xy=(band2[1], y_min), xytext=(0, 0), textcoords='offset points', fontsize=6)
    # else:
    # print(band1[0],band1[1])
    plt.plot([band1[0], band1[0]], [y_min * (1 + factor), y_max * (1 - factor)], color="blue")
    plt.plot([band1[1], band1[1]], [y_min * (1 + factor), y_max * (1 - factor)], color="blue")
    plt.annotate("VBM", xy=(band1[0], 0), xytext=(-12, 0), textcoords='offset points', fontsize=12, rotation=90)
    plt.annotate("CBM", xy=(band1[1], 0), xytext=(6, 0), textcoords='offset points', fontsize=12, rotation=90)
    print(fermi_level - vbm, temperature)

    if fermi_where:
        fermi_level = fermi_where + vbm
    plt.plot([fermi_level, fermi_level], [y_min, y_max], "b--")
    # label to the fermi_level
    # plt.annotate(f"{fermi_level:.2f} eV ({temperature} K)",
    #              xy=(fermi_level, y_min), xytext=(0, 0), textcoords='offset points', fontsize=8)

    plt.xlabel('Fermi level (eV)')
    plt.xlim(band1[0] - 0.1 * (band1[1] - band1[0]), band1[1] + 0.4 * (band1[1] - band1[0]))
    # plt.ylim(-2, 6)
    plt.ylabel('Energy (eV)')
    plt.legend(loc="upper right")
    plt.xticks(np.linspace(0, 5, 6) + band1[0], np.linspace(0, 5, 6))
    # plt.title("charge transition level of "+path)
    plt.savefig(path + ".png", dpi=300)
    plt.close('all')
    # plt.show()
    return None


if __name__ == '__main__':
    # a = {'vcs': [[-1, 1.823020758586241],
    #              [0, 1.3729962939999716],
    #              [-2, 5.404694152345081]],
    #      'vna': [[-1, 1.3735700405862288],
    #              [0, 0.8788255759999517],
    #              [-2, 4.950773434345122]],
    #      'vbi': [[-3, 7.281041466276626],
    #              [-2, 6.135769143345065],
    #              [-4, 11.125532718380692]],
    #      'vcl': [[1, -0.6121421954136753],
    #              [2, -0.6390288016548442],
    #              [0, 1.8334933399999596]],
    #      'bics': [[2, -1.045857132654921],
    #               [3, -0.8750348097233129],
    #               [1, 0.812154657586291],
    #               [0, 2.645685009000024]],
    #      'nacs': [[0, 0.9868907179999828],
    #               [1, 0.05460328458624419],
    #               [-1, 4.398005182586221]],
    #      'csna': [[0, 1.9567492819999766],
    #               [1, 1.7230237465862843],
    #               [-1, 5.358463746586235]],
    #      'bina': [[2, -1.8041478506549236],
    #               [3, -1.701345527723386],
    #               [1, 0.26802875558630623],
    #               [0, 2.616834291000024]],
    #      'csbi': [[-2, 6.119372849345124],
    #               [-1, 5.40484945558622],
    #               [-3, 9.799955172276647]],
    #      'nabi': [[-2, 4.221733567345133],
    #               [-1, 3.6161901735862263],
    #               [-3, 7.906865890276627]]}
    a = {"Al\-(Sc)": [[0, -4.8284], [-1, 7.469409962]],
         "Al\-(Ti)": [[1, -14.72332504], [0, -5.9061], [-1, 5.000104962]],
         "Al\-(V )": [[1, -14.55407504], [0, -7.118], [-1, 2.623754962]],
         "Al\-(Cr)": [[1, -14.33746504], [0, -8.2164], [-1, 2.289544962]],
         "Al\-(Mn)": [[1, -14.22778004], [0, -6.9167], [-1, 1.500559962]],
         "Al\-(Fe)": [[1, -10.72378004], [0, -5.1841], [-1, 3.682259962]],
         "Al\-(Co)": [[1, -7.264220038], [0, -1.5529], [-1, 6.237199962], [-2, 17.63687985]],
         "Al\-(Ni)": [[1, -3.375690038], [0, 2.8291], [-1, 9.631269962], [-2, 21.30171985]],
         "Al\-(Cu)": [[1, 2.968409962], [0, 7.9983], [-1, 15.14656996], [-2, 24.80131985]],
         "Al\-(Zn)": [[1, 4.865129962], [0, 10.161], [-1, 15.86114996]], }

    plot_ctl(a, 500, "Al2O3", [3.9873, 13.1415, 100])
