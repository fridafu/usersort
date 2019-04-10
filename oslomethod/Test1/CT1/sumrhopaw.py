import numpy as np
data = np.loadtxt('rhopaw.cnt')
data2 = np.loadtxt('rholev.cnt')
s = 0
s2 = 0
lvls = 0
# 28 is at
# 26 is at 3120 keV
# 20 is at 2400 keV
for i in range(0,30):
	s += data[i]
	s2 += data2[i]
	if data2[i] != 0:
		lvls += 1

print('rhopaw: ')
print(0.12*s)
print('rholev: ')
print(0.12*s2)
print('levels: ')
print(lvls)