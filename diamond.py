from pyphotonics.photoluminescence2 import Photoluminescence
import numpy as np
import matplotlib.pyplot as plt

if False:
    m=[78.971, 132.905, 118.711, 35.453]
    n=[1,8,3,24]
    mn=[]
    for i,j in zip(m,n):
        mn+=j*[i* 1.660539040e-27]

    path = "./"
    path_phonopy = "./"
    p = Photoluminescence("band.yaml","P1","P2",108, "phonopy", mn, 1000, shift_vector=[0, 0, 0])

else:
    m = np.zeros(63)
    for i in range(63):
        m[i] = 12.011 * 1.660539040e-27
    m[62] = 14.007 * 1.660539040e-27

    path = "./"
    path_phonopy = "./"
    #p = Photoluminescence("band1.yaml",path + "CONTCAR_GS",path + "CONTCAR_ES",189, "phonopy", m, 1000, shift_vector=[0.0, 0, 0])
    p = Photoluminescence("band1.yaml", path + "CONTCAR_GS", path + "CONTCAR_ES",  m)


print("Delta_R=", p.Delta_R)
print("Delta_Q=", p.Delta_Q)
print("HuangRhyes=", p.HuangRhyes)

plt.figure(figsize=(10, 10))
plt.plot(p.S_omega)
plt.ylabel('$S(\hbar\omega)$')
plt.xlabel('Phonon energy (meV)')
plt.xlim(0, 200)
# plt.ylim(0, 0.01)
plt.savefig('S_omega', bbox_inches='tight')

p.write_S('S')

A, I = p.PL(2, 2, 1.95)

plt.figure(figsize=(10, 10))
plt.plot(I.__abs__())
plt.ylabel('$I(\hbar\omega)$')
plt.xlabel('Photon energy (eV)')
#plt.xlim(1200, 2000)
x_values, labels = plt.xticks()
labels = [float(x)/p.resolution for x in x_values]
plt.xticks(x_values, labels)
#plt.ylim(0, 600)
plt.savefig('I', bbox_inches='tight')

plt.figure(figsize=(10, 10))
plt.plot(p.frequencies,p.S)
plt.savefig('S', bbox_inches='tight')