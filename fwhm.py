import numpy as np
import matplotlib.pyplot as plt

def fwhm(hw, s):
    def fg(t):
        return 2.36 * np.power(s, 0.5) * hw * np.power(np.tanh(hw / (2 * 8.6173303e-5 * t)), -0.5)
    return fg

a=[["SbCl$_6^{3-}$ ",0.723,31.7]   ,
["SbCl$_5^{2-}$ ",0.902,28.0]   ,
["BiCl$_6^{3-}$ cz ",0.099,31.5],
["BiCl$_6^{3-}$ ez ",0.273,31.5],
["BiCl$_5^{2-}$ ",0.481,29.0]   ,
["SeCl$_6^{2-}$ ",0.797,31.2]   ,
["TeCl$_6^{2-}$ cz ",0.234,31.1],
["TeCl$_6^{2-}$ ez ",0.575,31.1],]
a.sort(key=lambda x:x[1],reverse=True)
x=np.linspace(1,300,20)
for i in a:
    b=fwhm(i[2] / 1000, i[1]/i[2]*1000)
    plt.plot(x,b(x),label=i[0])
    #plt.annotate(i[0],xy=(x[-1],b(x[-1])))
plt.xlabel("Temperature (K)",fontsize=16)
plt.xlim(0,)
plt.ylabel("FWHM (eV)",fontsize=16)
plt.legend(fontsize=12)
#plt.show()
plt.savefig("fwhm.png",dpi=300)
plt.close(all)