import numpy as np
data = np.loadtxt('rhopaw.cnt')
s = 0

for i in range(1,30):
	s += data[i]
s = 1.2*s

print(s)