import MyModules as my
import numpy as np

def f(x):
    y1 = x[0][0] + x[1][0] + x[2][0]**3
    y2 = x[1][0]**2 + x[2][0]
    y3 = x[2][0]**3 + x[1][0]**2 + x[0][0]**3
    return np.array([y1, y2, y3]).reshape(3,1)
    
x = np.array([1, 2, 3]).reshape(3,1)

J = my.jacobian(f, x)