import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import linalg as LA
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib.widgets import Button

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
        step = 0.000001
        M = M + step * M_dot
        a = a + 1
    return result

limit = 0.01
M = [0.006, 0, 0]
B = [0, 0, 1]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d(-limit, limit)
ax.set_ylim3d(-limit, limit)
ax.set_zlim3d(-limit, limit) 

#code to create buttons
def Mup(event):
    global M
    M = M *1.1
    
def Mdown(event):
    global M
    M = M * 0.9
    
axprev = plt.axes([0.7, 0.85, 0.1, 0.075])
axnext = plt.axes([0.81, 0.85, 0.1, 0.075])
bnext = Button(axnext, 'M up')
bnext.on_clicked(Mup)
bprev = Button(axprev, 'M down')
bprev.on_clicked(Mdown)

def updateAnimation(M_init, B, t):
    ax.cla()
    ax.set_xlim3d(-limit, limit)
    ax.set_ylim3d(-limit, limit)
    ax.set_zlim3d(-limit, limit) 
    data = M_no_relax(M_init, B, t)
    global M 
    M = data[int((len(data)-1)*0.25)]
    data2 = np.transpose(np.reshape(data, (len(data),3), order='F'))
    return ax.plot(data2[0], data2[1], data2[2], label='Magnetization')

def init():
    return updateAnimation(M, B, 2000),

def animate(i):
    return updateAnimation(M, B, 2000),

updateAnimation(M, B, 2000)    

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=1, interval=2, blit=False)

plt.show()
    
    

