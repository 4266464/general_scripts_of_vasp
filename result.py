# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 17:17:13 2021

@author: dugue
"""
import numpy as np
import random
a=[["1Cs2AgInCl6_Bi",3.40,3.23,2.10,1.49],
["1Cs2KInCl6_Sb",3.87,3.85,2.50,3.00],
["1Cs2NaBiCl6_Sb",3.44,3.55,2.58,2.47],
["1Cs2NaInCl6_Sb",3.70,3.73,2.80,2.74],
["1Cs2NaLaCl6_Bi",3.76,3.72,3.63,3.58],
["1Cs2NaLaCl6_Sb",3.80,3.62,2.50,2.40],
["1Cs2NaScCl6_Sb",3.75,3.72,2.80,2.81], 
["1Cs2NaYCl6_Bi",3.80,3.78,3.72,3.68],
["1Cs2NaYCl6_Sb",3.75,3.67,2.70,2.62],
["1Cs2SnCl6_Bi+$V_{Cl}$",3.40,3.32,2.70,2.36],
["1KCl_Bi",3.70,3.70,3.20,3.21],
["1KCl_Sb",3.20,3.85,2.17,1.83],
["1Rb3InCl6_Sb",3.63,3.65,2.38,2.13]]

b=[["1Cs2AgInCl6_Bi",3.40,3.23,2.10,1.49,2.57351],
["1Cs2KInCl6_Sb",3.87,3.85,2.50,3.00,2.54712],
["1Cs2NaBiCl6_Sb",3.44,3.55,2.58,2.47,2.72536],
["1Cs2NaInCl6_Sb",3.70,3.73,2.80,2.74,2.55322],
["1Cs2NaLaCl6_Bi",3.76,3.72,3.63,3.58,2.77656],
["1Cs2NaLaCl6_Sb",3.80,3.62,2.50,2.40,2.77656],
["1Cs2NaScCl6_Sb",3.75,3.72,2.80,2.81,2.50912],
["1Cs2NaYCl6_Bi",3.80,3.78,3.72,3.68,2.64217],
["1Cs2NaYCl6_Sb",3.75,3.67,2.70,2.62,2.64217],
["1KCl_Bi",3.70,3.70,3.20,3.21,3.19180],
["1KCl_Sb",3.20,3.85,2.17,1.83,3.19180],
["1Rb3InCl6_Sb",3.63,3.65,2.38,2.13,2.54941]]

sb=[
["1Cs2NaBiCl6_Sb",3.44,3.55,2.58,2.47,2.72536],
["1Cs2NaInCl6_Sb",3.70,3.73,2.80,2.74,2.55322],
["1Cs2NaLaCl6_Sb",3.80,3.62,2.50,2.40,2.77656],
["1Cs2NaScCl6_Sb",3.75,3.72,2.80,2.81,2.50912],
["1Cs2NaYCl6_Sb",3.75,3.67,2.70,2.62,2.64217],
["1KCl_Sb",3.20,3.85,2.17,1.83,3.19180],
]

bi=[
["1Cs2NaLaCl6_Bi",3.76,3.72,3.63,3.58,2.77656],
["1Cs2NaYCl6_Bi",3.80,3.78,3.72,3.68,2.64217],
["1KCl_Bi",3.70,3.70,3.20,3.21,3.19180],
]

word_size=12

def scatter_plot_em(list_scat,n):
    import matplotlib.pyplot as plt
    for i in sorted(list_scat,key=lambda x:x[4]):
        if "_Bi" in i[0]:
            plt.scatter(i[4],i[3],color="red")
        else:
            plt.scatter(i[4],i[3],marker="^",color="blue")
        
        plt.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")[1:],
                     xy=(i[4],i[3]),xytext=(0,0),textcoords='offset points',fontsize=word_size)
    plt.annotate("Emission",
                     xy=(1.45,3.8),
                     xytext=(0,0),textcoords='offset points',fontsize=18)
            
        
    plt.plot([1.4,4],[1.4,4],linestyle="dashed",color="lightgray")
    plt.xlabel('$E$(calc) (eV)',fontsize=18)
    plt.xlim(1.4,4)
    plt.ylim(1.4,4)
    plt.ylabel('$E$(expt) (eV)',fontsize=18)
    plt.show()
    #plt.savefig("em"+str(n),dpi=300)
    plt.close("all")
    return None

def scatter_plot_ex(list_scat,n):
    import matplotlib.pyplot as plt
    prev=[]
    for i in sorted(list_scat,key=lambda x:x[2]):
        if "_Bi" in i[0]:
            plt.scatter(i[2],i[1],color="red")
        else:
            plt.scatter(i[2],i[1],marker="^",color="blue")
        
        if prev:
            if abs(prev[1]-i[1])<0.1 and abs(prev[0]-i[2])<0.07:
                #i[1]-=0.2*np.sign(prev[1]-i[1])
                #plt.plot([i[2],i[2]],[i[1],i[1]-0.5])
                d=0.4*random.random()
                if random.random()>0.5:
                    plt.plot([i[2],i[2]],[i[1],i[1]-d+0.03],"k--")
                    plt.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")[1:],
                     xy=(i[2],i[1]-d),xytext=(-10,0),textcoords='offset points',fontsize=word_size)
                else:
                    plt.plot([i[2],i[2]],[i[1],i[1]+d],"k--")
                    plt.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")[1:],
                     xy=(i[2],i[1]+d),xytext=(-10,0),textcoords='offset points',fontsize=word_size)
                    
            else:
                plt.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")[1:],
                     xy=(i[2],i[1]),xytext=(-10,10),textcoords='offset points',fontsize=word_size)
                prev=[i[2],i[1]]
        else:
            plt.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")[1:],
                 xy=(i[2],i[1]),xytext=(-10,10),textcoords='offset points',fontsize=word_size)
            prev=[i[2],i[1]]
    plt.annotate("Excitation",
                     xy=(3.15,3.8),
                     xytext=(0,0),textcoords='offset points',fontsize=18)
            
        
    plt.plot([3.1,4],[3.1,4],"k--")
    plt.xlabel('$E$(calc) (eV)',fontsize=18)
    plt.xlim(3.1,4)
    plt.ylim(3.1,4)
    plt.ylabel('$E$(expt) (eV)',fontsize=18)
    #plt.show()
    plt.savefig("ex"+str(n),dpi=300)
    plt.close("all")
    return None

# def scatter_plot_ex(list_scat):
#     import matplotlib.pyplot as plt
    
#     for i in list_scat:
#         flag=True
#         if "_Bi" in i[0]:
#             plt.scatter(i[2],i[1],color="red")
#         else:
#             plt.scatter(i[2],i[1],marker="^",color="blue")
#         for j in range(list_scat.index(i)):
#             if abs(i[1]-list_scat[j][1])<0.01 and abs(i[2]-list_scat[j][2])<0.3:
#                 flag=False
#         plt.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")[1:],xy=(i[2],i[1]),xytext=(0,0) if flag else (0,5),textcoords='offset points',fontsize=word_size)
            
        
#     plt.plot([3.1,4],[3.1,4],"k--")
#     plt.xlabel('Calculated $\Delta$Eexc (eV)')
#     plt.xlim(3.1,4)
#     plt.ylim(3.1,4)
#     plt.ylabel('Measured $\Delta$Eexc (eV)')
#     #plt.show()
#     plt.savefig("ex",dpi=300)
#     plt.close("all")
#     return None

def scatter_plot2(list_scat,n):
    import matplotlib.pyplot as plt
    import random
    import numpy as np
    label=["exp_ex","ex","exp_em","em"]
    for j in range(5)[-1:]:
        for i in list_scat:
            if "_Bi" in i[0]:
                plt.scatter(i[5],i[j],color="red")
                plt.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")[1:],xy=(i[5],i[j]),xytext=(-10,5),textcoords='offset points',fontsize=word_size)
            else:
                plt.scatter(i[5],i[j],marker="^",color="blue")
                plt.annotate(i[0].replace("Cs2","").replace("Cl6","").replace("_Sb","").replace("_Bi","")[1:],xy=(i[5],i[j]),xytext=(-10,5 if random.random()>0.5 else -15),textcoords='offset points',fontsize=word_size)
        x=np.linspace(2.4,3.3,2)
        coef=np.polyfit(list(zip(*sb))[5],list(zip(*sb))[4],1)
        plt.plot(x,np.poly1d(coef)(x),"k--")
        
        coef=np.polyfit(list(zip(*bi))[5],list(zip(*bi))[4],1)
        plt.plot(x,np.poly1d(coef)(x),"k--") 
        
        plt.xlabel('Bond length of $M$-Cl in host ($\AA$)',fontsize=16)
        plt.ylabel("$E_{emission}$(calc) (eV)",fontsize=16)
        
        # plt.xlabel('Structure',fontsize=10)
        # plt.ylabel("Emission",fontsize=14)
        # plt.xticks([])
        # plt.yticks([])
        
        plt.tight_layout()
        plt.show()
        #plt.savefig(label[j-1]+"priSb"+str(n),dpi=300)
        plt.close("all")
    return None

#scatter_plot2(b,1)
#scatter_plot_ex(a)

# for i in range(30):
#     #scatter_plot2(b,i)
scatter_plot_em(a,1)
    
    