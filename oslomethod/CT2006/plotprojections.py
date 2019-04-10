import numpy as np
from matplotlib import pyplot as plt
raw = np.loadtxt('projection_1072keV_raw.txt',dtype='float')
raw2 = raw[370:395]
resp = np.loadtxt('projection_1072keV_respmat.txt', dtype='float')
resp2 = resp[370:395]

n = 395-370#len(raw)
n = len(raw)
x = np.linspace(0,n,n)
xgs = 8*x - 2000 #correct for gain/shift
sr = np.sum(raw2/sum(raw))
sresp = np.sum(resp2/sum(resp))
print(sr)
print(sresp)
plt.plot(xgs,raw/sum(raw))
plt.plot(xgs,resp)
plt.xlabel(r'$\gamma$-energy')
plt.ylabel('Counts (normalized)')
#plt.xlim(960,1150)
plt.show()