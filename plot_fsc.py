# import numpy
# import Image
# import matplotlib.pyplot as plt

# bug =Image.open('C:/Users/sdh/Desktop/snow.jpg')
# arr =numpy.array(bug.getdata(),numpy.uint8).reshape(bug.size[1],bug.size[0],3)
# # plt.gray()
# plt.imshow(arr)
# plt.colorbar()
# plt.show()
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import random

# normal distribution center at x=0 and y=5
x=random.randint(1,100)
y=random.randint(1,100)
plt.plot(x, y)
# plt.colorbar()
plt.show()