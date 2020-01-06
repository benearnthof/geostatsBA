#!/usr/bin/env python
from __future__ import division
import numpy as np
from numpy import *
import matplotlib
#matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.special import jn, jn_zeros

def k(m,n):
    return jn_zeros(n,m)[m-1] # m is 0-indexed here

xs = linspace(-1, 1, 100)
ys = linspace(-1, 1, 100)
X, Y = meshgrid(xs, ys)
fns = 'cc'
f1 = {'s':sin,'c':cos}[fns[0]]
f2 = {'s':sin,'c':cos}[fns[1]]
n, m = 1, 1
v = 1
NT = 400
dt = 10/NT

def getZ(X, Y, i):
    global m,n,f1,f2,dt
    t = i * dt
    theta = arctan2(Y,X) # This does arctan(Y/X) but gets the sign right.
    R = sqrt(X**2 + Y**2)
    # We know z = J_n(k*r)*cos(n*theta)*cos(k*v*t)
    # 
    result = jn(n,k(m,n)*R)*f1(n*theta)*f2(k(m,n)*v*t)
    result[R>1] = 0  # we plot points from the square, but physically require this.
    return result


fig = plt.figure()
ax = plt.gca(projection='3d')#Axes3D(fig)
#wireframe = ax.plot_surface(X, Y, getZ(X,Y,0.0), rstride=2, cstride=2,cmap=cm.jet)
wireframe = ax.plot_surface(X, Y, getZ(X,Y,0.0), rstride=4, cstride=4, cmap=cm.jet, alpha=0.3)

ax.set_zlim(-1,1)
ax.axis("off")
ax.set_frame_on(False)

# precalculate Z

Z = [getZ(X,Y,i) for i in range(NT)]

def animate(i, ax, fig):
    global X,Y,Z
    global m,n
    #Z = getZ(X,Y,i)
    ax.cla()
    wireframe = ax.plot_surface(X, Y, Z[i], rstride=4, cstride=4, cmap=cm.jet, alpha=0.3)
    if n == 0:
        ax.set_zlim(-1,1)
    else:
        ax.set_zlim(-0.5,0.5)
    # ax = plt.gca(projection='3d'); ax._axis3don = False
    ax.axis("off")
    ax.set_frame_on(False)
    # frame.axes.get_xaxis().set_visible(False)
    # frame.axes.get_yaxis().set_visible(False)
    # ax.set_frame_on(False)
    #ax.view_init(elev=25., azim=140)
    return wireframe,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames=NT, fargs=(ax, fig), interval=2)
 
# this is how you save your animation to file:
anim.save('2d_sor_unstable.gif', writer='imagemagick', fps=30)
 
plt.show()

############## pseudo3d amination
#!/usr/bin/env python
from __future__ import print_function
"""
Based on the wireframe example script

If you make the plots, you can make the movie with

ffmpeg -q:a 5 -r 5 -b:v 19200 -i img%04d.png movie.mp4

or something similar.
"""
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np
from numpy import sin, cos, arctan, arctan2, array, sqrt, pi, linspace, meshgrid
import scipy
import scipy.special
from scipy.special import jn, jn_zeros
import time, sys, os

makeplots = True


v = 1
fns = 'cc'
ms = (1,2)
#ms = (1,2,3)
ns = (0,1)
#ns = (0,1,2)
    
f1 = {'s':sin,'c':cos}[fns[0]]
f2 = {'s':sin,'c':cos}[fns[1]]

def k(m,n):
    return jn_zeros(n,m)[m-1] # m is 0-indexed here


def generate(X, Y, t, m, n, v):
    theta = arctan2(Y,X)
    R = sqrt(X**2 + Y**2)
    # We know z = J_n(k*r)*      cos(n*theta)*cos(k*v*t)
    result =      jn(n,k(m,n)*R)*f1(n*theta)* f2(k(m,n)*v*t)
    result[R>1] = 0  # we plot points from the square, but physically require this.
    return result

plt.ion()
fig = plt.figure()
# It seems to me that things look best if you have 6 inches in width
# for each m and 2 in in height for each n.

fig.set_size_inches([6*len(ms),2*len(ns)])
axs = {}
rows, cols = len(ns), 2*len(ms)

idx = 1
for m in ms:
    axs[m] = {}
    for n in ns:
        axs[m][n] = (fig.add_subplot(rows,cols,idx, projection='3d'),
                    fig.add_subplot(rows,cols,idx+1))
        idx += 2
if 0:
    axs[1] = {0: (fig.add_subplot(rows,cols,1, projection='3d'), fig.add_subplot(rows,cols,2)),
              1: (fig.add_subplot(rows,cols,3, projection='3d'), fig.add_subplot(rows,cols,4)),}
    axs[2] = {0: (fig.add_subplot(rows,cols,5, projection='3d'), fig.add_subplot(rows,cols,6)),
              1: (fig.add_subplot(rows,cols,7, projection='3d'), fig.add_subplot(rows,cols,8)),}

xs = linspace(-1, 1, 100)
ys = linspace(-1, 1, 100)
X, Y = meshgrid(xs, ys)
#Z = generate(X, Y, 0.0)

wframes = {}
contours = {}
for m in ms:
    wframes[m] = {}
    contours[m] = {}
    for n in ns:
        wframes[m][n] = None
        contours[m][n] = None

if makeplots:
    dname = fns[0]+fns[1] + '.'.join([str(i) for i in ms]) + '_' + '.'.join([str(i) for i in ns])
    if not os.path.exists(dname):
        os.makedirs(dname)
tstart = time.time()

first = True
frames_per = 50
periods = 2
# Note: unlike sin and cos, Jn's zeros are not integer multiples of
# each other.  Therefore, this loop goes over a defined number of
# periods of the lowest mode, but others won't fit evenly.
junk = input('go for barney')
for (idx,t) in enumerate(linspace(0, periods*2*pi/jn_zeros(0,1)[0], periods*frames_per)):
    if not first:
        for m in ms:
            for n in ns:
                axs[m][n][0].collections.remove(wframes[m][n])
                for c in contours[m][n].collections:
                    axs[m][n][1].collections.remove(c)
                    
    first = False

    for m in ms:
        for n in ns:
            Z = generate(X, Y, t, m, n, v)
            wframes[m][n] = axs[m][n][0].plot_surface(X, Y, Z, cmap=cm.coolwarm, rstride=4, cstride=4, alpha=0.3)

            if n == 0:
                axs[m][n][0].set_zlim(-0.7,.7)
                axs[m][n][1].imshow(Z,vmin=-0.7,vmax=0.7)
            else:
                axs[m][n][0].set_zlim(-0.5,0.5)
                axs[m][n][1].imshow(Z,vmin=-0.5,vmax=0.5)
            # The funny business with levels here is because you won't
            # get a contour exactly at zero that necessarily tracks
            # around both sides of the circle due to the fact that
            # we've discretized things.
            levels = [-0.000000001,0.0,0.000000001]
            contours[m][n] = axs[m][n][1].contour(Z, levels, colors='k',
                                                  linestyles='solid', linewidths=2)
    plt.draw()
    if makeplots:
        plt.savefig(dname + '/img' + '%04d'%(idx) + '.png')

print ('FPS: %f' % (periods*frames_per / (time.time() - tstart)))
