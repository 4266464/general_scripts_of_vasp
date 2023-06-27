# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 21:30:45 2021

@author: dugue
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 06:40:18 2021

@author: dugue
"""

a=[["Cs2AgInCl6_Bi",2.57351,2.67313,2.607476667,2.82458997,0],
["Cs2NaInCl6_Bi",2.55322,2.68561,2.57526,2.825738613,2.721544553],
["Cs2NaLaCl6_Bi",2.77656,2.713906667,2.591783333,2.877410074,2.764377423],
["Cs2NaScCl6_Bi",2.50912,2.67775,2.569533333,2.808101905,2.709274348],
["Cs2NaYCl6_Bi",2.64217,2.697876667,2.58193,2.848635431,2.739110451],
["Cs2SnCl6_Bi",2.47435,2.70157,2.587913333,2.88478409,2.758610203],
["Cs2ZrCl6_Bi",2.48548,2.704173333,2.591463333,2.887376379,2.768701052],
["KCl_Bi",3.19180,2.73688,2.60041,2.975288949,2.814015061],
["Cs2AgInCl6_Sb",2.57351,2.60262,2.488469869,2.748445123,0],
["Cs2KInCl6_c2h_Sb",2.54440,2.626722078,2.487359732,2.853403717,2.6674902],
["Cs2NaBiCl6_Sb",2.72536,2.64774,2.50091,2.801228204,2.653882285],
["Cs2NaInCl6_Sb",2.55322,2.62807,2.486276667,2.765786811,2.623862952],
["Cs2NaLaCl6_Sb",2.77656,2.648626667,2.501193333,2.813027554,2.657399602],
["Cs2NaScCl6_Sb",2.50912,2.622533333,2.48536,2.752592854,2.614542855],
["Cs2NaYCl6_Sb",2.64217,2.64052,2.494906667,2.785560894,2.638453908],
["Cs2SnCl6_Sb",2.47435,2.63494,2.49528,0,2.665019409],
["Cs2ZrCl6_Sb",2.48548,2.641936667,2.498696667,2.837962867,2.672769905],
["KCl_Sb",3.19180,2.66667,2.503836667,2.91455163,2.723908911],
["Rb3InCl6_c2h_Sb",2.54941,2.632385289,2.493248532,2.832466528,2.668789594],
["Cs2KInCl6_oh_Sb",2.54712,2.625733333,2.489783333,2.751773333,2.614539201]]

b=[["Cs2NaInCl6_Bi",2.55322,2.68561,2.57526,2.825738613,2.721544553],
["Cs2NaLaCl6_Bi",2.77656,2.713906667,2.591783333,2.877410074,2.764377423],
["Cs2NaScCl6_Bi",2.50912,2.67775,2.569533333,2.808101905,2.709274348],
["Cs2NaYCl6_Bi",2.64217,2.697876667,2.58193,2.848635431,2.739110451],
["KCl_Bi",3.19180,2.73688,2.60041,2.975288949,2.814015061],
["Cs2NaBiCl6_Sb",2.72536,2.64774,2.50091,2.801228204,2.653882285],
["Cs2NaInCl6_Sb",2.55322,2.62807,2.486276667,2.765786811,2.623862952],
["Cs2NaLaCl6_Sb",2.77656,2.648626667,2.501193333,2.813027554,2.657399602],
["Cs2NaScCl6_Sb",2.50912,2.622533333,2.48536,2.752592854,2.614542855],
["Cs2NaYCl6_Sb",2.64217,2.64052,2.494906667,2.785560894,2.638453908],
["KCl_Sb",3.19180,2.66667,2.503836667,2.91455163,2.723908911],
["Cs2KInCl6_oh_Sb",2.54712,2.625733333,2.489783333,2.751773333,2.614539201]]

sd={"Cs2KInCl6_oh_Sb":0.327,
"Cs2NaBiCl6_Sb":0.437,
"Cs2NaInCl6_Sb":0.377,
"Cs2NaLaCl6_Sb":0.472,
"Cs2NaScCl6_Sb":0.344,
"Cs2NaYCl6_Sb":0.403,
"Cs2SnCl6_Sb":0.504,
"Cs2ZrCl6_Sb":0.507,
"KCl_Sb":0.75,
"Rb3InCl6_oh_Sb":0.317}

sda={"Cs2KInCl6_oh_Sb":0.327,
"Cs2NaBiCl6_Sb":0.437,
"Cs2NaInCl6_Sb":0.377,
"Cs2NaLaCl6_Sb":0.472,
"Cs2NaScCl6_Sb":0.344,
"Cs2NaYCl6_Sb":0.403,
"KCl_Sb":0.75,
"Rb3InCl6_oh_Sb":0.317}

sa=[i for i in a if "_Sb" in i[0]]
ba=[i for i in a if "_Bi" in i[0]]



s=[i for i in b if "_Sb" in i[0]]
b=[i for i in b if "_Bi" in i[0]]


def scatter_plot_ex(list_scat):
    import matplotlib.pyplot as plt
    import numpy as np
    fig,ax=plt.subplots()
    for j,k in zip(range(6)[2:],["grd","m1","p1","sp"]):
        for i in list_scat:
            if i[j]:
                ax.scatter(i[1],i[j],color="red",s=1)
                ax.annotate(i[0].replace("Cs2","").replace("Cl6","").split("_")[0]+"_"+k,xy=(i[1],i[j]),xytext=(0,0),
                              textcoords='offset points',fontsize=5)
    
        x=np.linspace(2.4,3.3,2)
        coef=np.polyfit(list(zip(*s))[1],list(zip(*s))[j],1)
        ax.plot(x,np.poly1d(coef)(x),"k--")
        
        # coef=np.polyfit(list(zip(*b))[1],list(zip(*b))[j],1)
        # plt.plot(x,np.poly1d(coef)(x),"k--")
    ax2=ax.twinx()
    xsd=[]
    ysd=[]
    for i in list_scat:
        if i[0] in sd:
            ax2.scatter(i[1],sd[i[0]],color="red",s=1)
            if i[0] in sda:
                xsd.append(i[1])
                ysd.append(sd[i[0]])
            ax2.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")+"_di",xy=(i[1],sd[i[0]]),xytext=(0,0),
                         textcoords='offset points',fontsize=5)
    x=np.linspace(2.4,3.3,2)
    coef=np.polyfit(xsd,ysd,1)
    ax2.plot(x,np.poly1d(coef)(x),"k--")

    ax.plot([2.4,3.2],[2.4,3.2],linestyle=(0,(1,10)),color="black")
    #ax.xlabel('pri')
    #ax.xlim(2.4,3.2)
    #plt.ylim(2.4,3.2)
    #plt.ylabel(k)
    #plt.savefig("1001s",dpi=300)
    plt.show()
    plt.close("all")

    #plt.plot([2.4,3.2],[2.4,3.2],"k--")
    #plt.xlabel('pri')
    #plt.xlim(2.4,3.2)
    #plt.ylim(2.4,3.2)
    #plt.ylabel('grd')
    
    # x=np.linspace(2.4,3.2,2)
    # coef=np.polyfit(list(zip(*s))[1],list(zip(*s))[2],1)
    # plt.plot(x,np.poly1d(coef)(x),"k--")
    
    # coef=np.polyfit(list(zip(*b))[1],list(zip(*b))[2],1)
    # plt.plot(x,np.poly1d(coef)(x),"k--")
    
    # coef=np.polyfit(list(zip(*b2))[1],list(zip(*b2))[2],1)
    # plt.plot(x,np.poly1d(coef)(x),"k--")
    
    return None


scatter_plot_ex(sa)