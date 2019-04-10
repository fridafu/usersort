from matplotlib import pyplot as plt
import numpy as np
x = np.linspace(0,10,10000)
y = 0.02*8**x
plt.semilogy(x,y)
plt.show()