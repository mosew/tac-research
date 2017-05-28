import numpy as np
from scipy.linalg import expm

x = np.array([[1,1]])
x.shape

y = np.array([2,2])

z = np.dot(x,y)
print z