# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 23:46:31 2021

@author: dugue
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_ctl(band_range, charge_list, temp, path):
    x = np.linspace(*band_range, 100)
    fermi = []
    for i in charge_list:
        y = []
        for j in i:
            y.append([[j[0] * k + j[1], j[0]] for k in x])

        node_list = []
        y_min = []
        for j in range(len(x)):
            min_value = [1e5]
            for k in range(len(y)):
                if y[k][j][0] < min_value[0]:
                    min_value = y[k][j]
            if y_min and not min_value[1] == y_min[-1][1]:
                node_list.append((x[j], (min_value[0] + y_min[-1][0]) * 0.5))
            y_min.append(min_value)
        fermi.append(y_min)

        # y_min_max=list(zip(*y_min))[0]
        # y_min_min = min(y_min_max)
        # y_min_max = max(y_min_max)

        plt.plot(x, list(zip(*y_min))[0], label=i[1][-1])
        # plt.plot([band_range[0], band_range[1]], [y_min[0][0], y_min[-1][0]], "ro")
        plt.annotate(i[1][-1], xy=(x[0], y_min[0][0]), xytext=(0, 0), textcoords='offset points', fontsize=6)
        # plt.annotate(i[1][-1], xy=(x[-1], y_min[-1][0]),xytext=(0, 0), textcoords='offset points', fontsize=6)
        # plt.plot([band_range[0], band_range[0]], [y_min_min *(1+factor), y_min_max*(1-factor)], color="blue")
        # plt.plot([band_range[1], band_range[1]], [y_min_min *(1+factor), y_min_max*(1-factor)], color="blue")

        # for j in node_list:
        #     plt.plot(j[0], j[1], "ro")
        # plt.annotate(i[1][-1], xy=(j[0], j[1]), xytext=(0, 0),textcoords='offset points', fontsize=6)

    fermi = list(zip(*fermi))
    y_max = list(zip(*fermi[0]))[0] + list(zip(*fermi[-1]))[0]
    y_min = min(y_max)
    y_max = max(y_max)
    factor = 2e-4
    plt.plot([band_range[0], band_range[0]], [y_min * (1 + factor), y_max * (1 - factor)], color="blue")
    plt.plot([band_range[1], band_range[1]], [y_min * (1 + factor), y_max * (1 - factor)], color="blue")
    plt.annotate("VBM", xy=(x[0], y_min), xytext=(0, 0), textcoords='offset points', fontsize=6)
    plt.annotate("CBM", xy=(x[-1], y_min), xytext=(0, 0), textcoords='offset points', fontsize=6)
    for i, j in zip(fermi, x):
        # print(j,sum([k[1]*np.exp(-k[0]/temp) for k in i]))
        if sum([k[1] * np.exp(-k[0] / temp) for k in i]) < 0:
            plt.plot([j, j], [y_min, y_max], "b--")
            plt.annotate("Fermi level=" + str(round(j, 2)) + "eV under " + str(round(temp / 0.025852 * 300)) + "K",
                         xy=(j, y_min), xytext=(0, 0), textcoords='offset points', fontsize=6)
            break

    plt.xlabel('band(eV)')
    plt.xlim(band_range[0] - 0.1 * (band_range[1] - band_range[0]),
             band_range[1] + 0.3 * (band_range[1] - band_range[0]))
    # plt.ylim(y_min_min*(1+factor),)
    plt.ylabel('energy(eV)')
    plt.legend(loc="upper right")
    # plt.title("charge transition level of "+path)
    # plt.legend()
    # plt.show()
    plt.savefig(path + ".png", dpi=300)
    plt.close('all')
    return list(zip(x, y))


a = [[[-1, 1.823020758586241, 'vcs'],
      [0, 1.3729962939999716, 'vcs'],
      [-2, 5.404694152345081, 'vcs']],
     [[-1, 1.3735700405862288, 'vna'],
      [0, 0.8788255759999517, 'vna'],
      [-2, 4.950773434345122, 'vna']],
     [[-3, 7.281041466276626, 'vbi'],
      [-2, 6.135769143345065, 'vbi'],
      [-4, 11.125532718380692, 'vbi']],
     [[1, -0.6121421954136753, 'vcl'],
      [2, -0.6390288016548442, 'vcl'],
      [0, 1.8334933399999596, 'vcl']],
     [[2, -1.045857132654921, 'bics'],
      [3, -0.8750348097233129, 'bics'],
      [1, 0.812154657586291, 'bics'],
      [0, 2.645685009000024, 'bics']],
     [[0, 0.9868907179999828, 'nacs'],
      [1, 0.05460328458624419, 'nacs'],
      [-1, 4.398005182586221, 'nacs']],
     [[0, 1.9567492819999766, 'csna'],
      [1, 1.7230237465862843, 'csna'],
      [-1, 5.358463746586235, 'csna']],
     [[2, -1.8041478506549236, 'bina'],
      [3, -1.701345527723386, 'bina'],
      [1, 0.26802875558630623, 'bina'],
      [0, 2.616834291000024, 'bina']],
     [[-2, 6.119372849345124, 'csbi'],
      [-1, 5.40484945558622, 'csbi'],
      [-3, 9.799955172276647, 'csbi']],
     [[-2, 4.221733567345133, 'nabi'],
      [-1, 3.6161901735862263, 'nabi'],
      [-3, 7.906865890276627, 'nabi']]]

# b=[-1.04494086,4.23386505]
band_range = [0.26434953, 3.29911075]

plot_ctl(band_range, a, 0.05, "fh")
