import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from numpy import linalg as LA

# gamma / 2 Pi in n MHz·T−1
gamma = 42.576

# computes gamma dot without relaxation
def M_dot_no_relax(M,B):
    return gamma * np.cross(M,B) * math.pi * 2

def M_no_relax(M_init, B, t):
    result = []
    a = 0;
    M = M_init;
    while a < t:
        result.append(M)
        #print ("M " + str(M))
        M_dot = M_dot_no_relax(M, B)
        #print ("M_dot", M_dot, str(LA.norm(M - M_dot)))
        step = 0.001# * LA.norm(M - M_dot)
        M = M + step * M_dot
        a = a + 1
    return result

M = [0.006, 0, 0]
B = [0, 0, 1]
data = M_no_relax(M, B, 200)
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
data2 = np.transpose(np.reshape(data, (len(data),3), order='F'))
plt.plot(data2[0], data2[1])
plt.show()
    
