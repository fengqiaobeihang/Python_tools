import numpy as np
import matplotlib.pylab as plt
d = 0.1
ds = np.array([150,200,210,220,230])
L = 0.334
P = 0.02
HI = -L*(ds*d-P*d*1000)/(ds*d)
plt.plot(HI)
plt.show()
H = sum(HI)/5
print (H)
print (HI)