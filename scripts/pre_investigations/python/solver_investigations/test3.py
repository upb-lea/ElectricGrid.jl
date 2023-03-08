import numpy as np

a = np.random.randn(2000, 2000)
x = np.random.randn(2000)

for i in range(50000):
    ax = np.dot(a, x)