from numpy import *
from matlotlib import pyplot as plt

N = 10000
a_0 = 5.2e-11 # bohr radius m
e = 1.6e-19 #Coulombs
x = np.zeros(N,N)
v = np.zeros(N,N)
t = np.zeros(N)
m = 9.1e-31 # electron mass kg
dt = 1
F = (k*e*-e)/(r**2)#coulomb k*

#a = np.zeros(N,N)
x[0,0] = (-500*a_0, 100*a_0) 
for i in range(N):
	a = F/m
	x[i+1] = x[i] + v[i+1]*dt
	v[i+1] = v[i] + a*dt
	#a(i+1) = 
	t[i+1] = t[i] + dt

plot(x,t)