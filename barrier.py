# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 15:47:03 2021

@author: dugue
"""

import matplotlib.pyplot as plt
import numpy as np

a=[[0,-170.91222,-172.39887,"vb",],
[0.1,-170.88999,-172.61644,"",],
[0.2,-170.8238,-172.79052,"",],
[0.3,-170.71745,-172.92374,"",],
[0.4,-170.54416,-173.01594,"",],
[0.5,-170.38826,-173.06811,"",],
[0.6,-170.17501,-173.08181,"",],
[0.7,-170.53173,-173.05651,"",],
[0.8,-170.60556,-172.99229,"",],
[0.9,-170.66563,-172.88978,"",],
[1,-170.68337,-172.74761,"ste",],
[1.1,-170.68068,-172.81877,"",],
[1.2,-170.67456,-172.88152,"",],
[1.3,-170.66708,-172.93615,"",],
[1.4,-170.66013,-172.98251,"",],
[1.5,-170.65435,-173.01993,"",],
[1.6,-170.65149,-173.049,"",],
[1.7,-170.6526,-173.0694,"",],
[1.8,-170.65509,-173.08064,"",],
[1.9,-170.65768,-173.08246,"",],
[2,-170.65944,-173.07427,"cb",],
[2.1,-170.65553,-173.147,"",],
[2.2,-170.64399,-173.21221,"",],
[2.3,-170.62498,-173.27005,"",],
[2.4,-170.59932,-173.32121,"",],
[2.5,-170.56612,-173.3647,"",],
[2.6,-170.52599,-173.40104,"",],
[2.7,-170.47849,-173.42985,"",],
[2.8,-170.42308,-173.45037,"",],
[2.9,-170.36022,-173.46304,"",],
[3,-170.28956,-173.46751,"grd",],
[3.1,-170.32321,-173.4567,"",],
[3.2,-170.39598,-173.42317,"",],
[3.3,-170.46013,-173.36751,"",],
[3.4,-170.54124,-173.29007,"",],
[3.5,-170.65384,-173.19237,"",],
[3.6,-170.74631,-173.07376,"",],
[3.7,-170.81882,-172.93511,"",],
[3.8,-170.87096,-172.77644,"",],
[3.9,-170.90244,-172.59791,"",],
[4,-170.91222,-172.39887,"vb",]]



y=list(zip(*a))[2]
ymax=max(y)
ymin=min(y)

f,(ax,ax2)=plt.subplots(2,1,sharex=True)

for i in a:
    ax2.scatter(i[0],i[2]-ymin,color="red")
    ax.scatter(i[0],i[1]-ymin,color="blue")
        

ax.set_ylim(2.4,3.4)
ax2.set_ylim(-0.2,1.2)
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)
ax2.xaxis.tick_bottom()
d=.015
kwargs=dict(transform=ax.transAxes,color='k',clip_on=False)
ax.plot((-d,+d),(-d,+d),**kwargs)
ax.plot((1-d,1+d),(-d,+d),**kwargs)

kwargs.update(transform=ax2.transAxes)
ax2.plot((-d,+d),(1-d,1+d),**kwargs)
ax2.plot((1-d,1+d),(1-d,1+d),**kwargs)
#plt.subplots_adjust(left=0.125,bottom=0.1,right=0.9,top=0.9,wspace=0.2,hspace=0.35)
plt.subplots_adjust(hspace=0.1)
for i in range(5):
    ax.plot([i,i],[2.4,3.4],"k--")
    ax2.plot([i,i],[-0.2,1.2],"k--")
plt.xticks([0,1,2,3,4],["CT","STE","MMCT","grd","CT"])
plt.xlabel("Geometry configuration")
plt.ylabel("Energy (eV)")
#plt.show()
plt.savefig("barrier_agin.png",dpi=300)



# y=list(zip(*a))[1]
# ymax=max(y)
# ymin=min(y)

# # f,(ax,ax2)=plt.subplots(2,1,sharex=True)

# for i in a:
#     if len(i)==3:
#         if i[2]=="cb":
#             plt.scatter(i[0],i[1]-ymin,color="red")
#         else:
#             plt.scatter(i[0],i[1]-ymin,color="blue")
#     else:
#         plt.scatter(i[0],i[1]-ymin,color="green")
# # ax.set_ylim(.78, 1.)
# # ax2.set_ylim(.78, 1.)
# for i in [0,1.2537196171127465,2.16477932,3.03948926]:
#     plt.plot([i,i],[0,ymax-ymin],"b--")
# plt.xticks([0,1.2537196171127465,2.16477932,3.03948926],["cb","ez","grd","cb"])
# plt.savefig("zzz.png",dpi=300)

