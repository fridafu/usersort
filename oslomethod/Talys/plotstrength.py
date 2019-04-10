from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt('transext.nrm',dtype=float)

a0 =  -0.8360;
a1 =   0.1280;
energy = np.zeros(len(data))
trans = np.zeros(len(data))
for i in range(len(data)):
    energy[i] = a0 + a1*i;
    trans[i] = data[i]/(2.*3.14*energy[i]*energy[i]*energy[i])


plt.plot(energy,trans)
#plt.show()

data2 = np.loadtxt('UnfoldedCrossSections_evaluated_at_Emax_RSFZn68.dat',dtype = float)
plt.plot(data2[:,0], data2[:,1])
plt.show()