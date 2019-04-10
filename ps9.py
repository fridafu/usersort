from numpy import random
N_p = 100000 # number of photon packets
L_s = 10e43 #erg/s luminosity of supernova
L_p = L_s/N_p
tau =  0.01
c = 3e10 #cm/s
beta = np.sqrt((1-gamma**-2))
gamma = 10#1/sqrt(1-beta**2)
E_in = 1 #eV

for i in range(N_p):
	costheta_in = random.uniform(-1,1)
	costheta_out = random.uniform(-1,1)
	E_out = E_in * gamma**2 * (1-beta*costheta_in)*(1+beta*costheta_out)
	L_out = L_p*(E_out/E_in)*(1-np.exp(-tau))

#plt.loglog(E_out,L_out)
plt.hist(L_out,N_p)
plt.show()