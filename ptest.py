# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 22:30:12 2021

@author: dugue
"""

import matplotlib.pyplot as plt
from adjustText import adjust_text
from scipy import interpolate
import numpy as np

together = [(0, 1.0, 0.4), (25, 1.0127692669427917, 0.41), (50, 1.016404709797609, 0.41), (75, 1.1043426359673716, 0.42), (100, 1.1610446924342996, 0.44), (125, 1.1685687930691457, 0.43), (150, 1.3486407784550272, 0.45), (250, 1.4013999168008104, 0.45)]
together.sort()

text = [x for (x,y,z) in together]
eucs = [y for (x,y,z) in together]
covers = [z for (x,y,z) in together]
def plot_eucs_covers():
    plt.plot(eucs,covers,color="black", alpha=0.5)
    texts = []
    for xt, yt, s in zip(eucs, covers, text):
        texts.append(plt.text(xt, yt, s))
    return texts

texts = plot_eucs_covers()
f = interpolate.interp1d(eucs, covers)
x = np.linspace(min(eucs), max(eucs), 500)
y = f(x)
#adjust_text(texts, x, y, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
adjust_text(texts, x, y)
plt.text(1.2,0.41,"lalala")