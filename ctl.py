# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 05:02:47 2021

@author: dugue
"""

import matplotlib.pyplot as plt
# ctl np
# a=[["1Cs2NaLaCl6_Bi",3.629340507,3.65732],
# ["1Cs2NaLaCl6_Sb",3.108160507,3.09592],
# ["1KCl_Sb",2.869782772,2.83038],
# ["1Cs2NaBiCl6_Sb",3.120912519,3.07402],
# ["1Rb3InCl6oh_Sb",3.685885231,3.61746],
# ["1KCl_Bi",3.667232772,3.58947],
# ["1Cs2NaYCl6_Bi",3.857671648,3.73925],
# ["1Cs2ZrCl6_Sb",3.117571858,2.99651],
# ["1Cs2NaYCl6_Sb",3.369681648,3.21237],
# ["1Cs2ZrCl6_Bi",3.756071858,3.5706],
# ["1Cs2KInCl6oh_Sb",3.723142137,3.51547],
# ["1Cs2HfCl6_Bi",3.796051044,3.58306],
# ["1Cs2NaInCl6_Bi",4.047605452,3.8226],
# ["1Cs2NaScCl6_Bi",4.091731473,3.83283],
# ["1Cs2SnCl6_Bi",3.876404836,3.6122],
# ["1Cs2NaInCl6_Sb",3.586105452,3.3177],
# ["1Cs2NaScCl6_Sb",3.641941473,3.3462],
# ["1Cs2ZrCl6_vcl_Bi",3.244571858,2.79539],
# ["1Cs2HfCl6_vcl_Bi",3.274421044,2.80259],
# ["1Cs2ZrCl6_vcl_Sb",2.835241858,2.30343],
# ["1Cs2SnCl6_vcl_Bi",3.364504836,2.82485],
# ["1Cs2SnCl6_vcl_Sb",2.968634836,2.33816]]

a=[["1Cs2NaLaCl6_Bi",3.629340507,3.65732],
["1Cs2NaLaCl6_Sb",3.108160507,3.09592],
["1KCl_Sb",2.869782772,2.83038],
["1Cs2NaBiCl6_Sb",3.120912519,3.07402],
["1Rb3InCl6oh_Sb",3.685885231,3.61746],
["1KCl_Bi",3.667232772,3.58947],
["1Cs2NaYCl6_Bi",3.857671648,3.73925],
["1Cs2NaYCl6_Sb",3.369681648,3.21237],
["1Cs2KInCl6oh_Sb",3.723142137,3.51547],
["1Cs2HfCl6_Bi",3.796051044,3.58306],
["1Cs2NaInCl6_Bi",4.047605452,3.8226],
["1Cs2NaScCl6_Bi",4.091731473,3.83283],
["1Cs2SnCl6_Bi",3.876404836,3.6122],
["1Cs2NaInCl6_Sb",3.586105452,3.3177],
["1Cs2NaScCl6_Sb",3.641941473,3.3462],
["1Cs2SnCl6_vcl_Bi",3.364504836,2.82485],
["1Cs2SnCl6_vcl_Sb",2.968634836,2.33816]]

s=[["1Cs2NaBiCl6_Sb",3.120912519,3.07402],
["1Cs2NaInCl6_Sb",3.586105452,3.3177],
["1Cs2NaLaCl6_Sb",3.108160507,3.09592],
["1Cs2NaScCl6_Sb",3.641941473,3.3462],
["1Cs2NaYCl6_Sb",3.369681648,3.21237]]

b=[["1Cs2NaInCl6_Bi",4.047605452,3.8226],
["1Cs2NaLaCl6_Bi",3.629340507,3.65732],
["1Cs2NaScCl6_Bi",4.091731473,3.83283],
["1Cs2NaYCl6_Bi",3.857671648,3.73925]]

b2=[["1Cs2HfCl6_Bi",3.796051044,3.58306],
["1Cs2SnCl6_Bi",3.876404836,3.6122],
["1Cs2ZrCl6_Bi",3.756071858,3.5706]]



def scatter_plot_ex(list_scat):
    import matplotlib.pyplot as plt
    import numpy as np
    
    for i in list_scat:
        flag=True
        if "_Bi" in i[0]:
            plt.scatter(i[1],i[2],color="red")
        else:
            plt.scatter(i[1],i[2],marker="^",color="blue")
        for j in range(list_scat.index(i)):
            if abs(i[1]-list_scat[j][1])<0.3 and abs(i[2]-list_scat[j][2])<0.01:
                flag=False
        plt.annotate(i[0][1:].replace("Cs2","").replace("Cl6","").replace("_Bi","").replace("_Sb",""),xy=(i[1],i[2]),xytext=(0,0) if flag else (0,5),textcoords='offset points',fontsize=10)
            
        
    plt.plot([2.7,4.2],[2.7,4.2],"k--")
    plt.xlabel('$iv_{ZPL}$ (eV)')
    plt.xlim(2.7,4.2)
    plt.ylim(2.7,4.2)
    plt.ylabel('ZPL (eV)')
    
    x=np.linspace(2.7,4.2,2)
    coef=np.polyfit(list(zip(*s))[1],list(zip(*s))[2],1)
    plt.plot(x,np.poly1d(coef)(x),"k--")
    
    coef=np.polyfit(list(zip(*b))[1],list(zip(*b))[2],1)
    plt.plot(x,np.poly1d(coef)(x),"k--")
    
    # coef=np.polyfit(list(zip(*b2))[1],list(zip(*b2))[2],1)
    # plt.plot(x,np.poly1d(coef)(x),"k--")
    
    
    #plt.savefig("1001",dpi=300)
    plt.show()
    plt.close("all")
    return None

scatter_plot_ex(a)

a.sort(key=lambda x:x[1])
for i in a:
    if "_Sb" in i[0]:
        print(i)