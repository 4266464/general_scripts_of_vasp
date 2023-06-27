import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = plt.axes(projection="3d")
peak = [[0, 0, 0], [1, 1.73205, 0], [-1, 1.73205, 0], [0, 1.1547, 1.63299]]
big = [[0, 0.515260243, 0.728686953],
       [0.553771332, 1.474419879, 0.728686953],
       [-0.553771332, 1.474419879, 0.728686953]]
zero = [[0, 0, 0], [1, 1.73205, 0], [-1, 1.73205, 0]]
#    [[0, 0.429256506, 0.60705948],
# [0.628252788, 1.517421747, 0.60705948],
# [-0.628252788, 1.517421747, 0.60705948]]
small = [[0, 1.051678439, 1.487295725],
         [0.089219331, 1.206210781, 1.487295725],
         [-0.089219331, 1.206210781, 1.487295725]]
line = [[0, 0.866025, 0], [0.008921933, 1.206210781, 1.487298605]]
for i in [peak, big, small]:
    for j in range(len(i)):
        for k in range(len(i))[j + 1:]:
            ax.plot3D(*list(zip(i[j], i[k])), "k")
for i in [line]:
    for j in range(len(i)):
        for k in range(len(i))[j + 1:]:
            ax.plot3D(*list(zip(i[j], i[k])), "cyan", linewidth=4)

for i in [big, small, zero]:
    ax.plot_trisurf(*list(zip(i[0], i[1], i[2])), alpha=0.5, color="orangered")

label = {'$\mu_{Al}$': [0, 1.932050808, 0],
         '$\mu_{Mg}$': [0.673205, 0.6, 0],
         '$\mu_{Si}$': [-0.673205, 0.766025, 0],
         '$\mu_O$': [0, 0.42735, 0.966495]}
for i in label:
    # ax.scatter(*label[i],color="k")
    ax.text(*label[i], i, fontsize=20)

label = {'-4.9 eV': [0, 0.5, 1.4],
         '-2.4 eV': [0, 0, 0.6],
         '0.0 eV': [0, -0.3, 0]}
for i in label:
    # ax.scatter(*label[i],color="k")
    ax.text(*label[i], i, fontsize=16)

ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
plt.axis("off")
plt.show()

# para=[[0.120299, 0.0505137, 0, -0.1, 0, 0],
#  [0.120314, 0.0836698, 0.31398, 0.3, 1.376095226, 0],
#  [0.193998, 0.041308, 0.40491, 0.5, 0.58421226, 1.821236482]]
# para = [[1, 0, 0], [-1, 0, 0], [0, 1.732050808, 0]]
# length=50
# x = np.linspace(-4, 4, length)
# y = np.linspace(-4, 4, length)
# x, y = np.meshgrid(x, y)
# fig = plt.figure()
# ax = plt.axes(projection="3d")
# for i in para:
#     z = ((x - i[0])) ** 2  + ((y - i[1])) ** 2+i[2]
#     #if ((x - i[0])) ** 2  + ((y - i[1])) ** 2+i[2]<8 else 1e8
#     for k in range(20):
#         for j in range(20):
#             if z[k][j]>8:
#                 z[k][j]=0
#     print(z)
#
#     ax.plot_surface(x, y, z, alpha=0.5)
# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Z")
# ax.set_zlim(0,8)
# plt.show()
