from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import matplotlib.pyplot as plt
import numpy as np

W_FILE_PATH = '/home/glen/Projects/ice_deflection/deflections.txt'

_w_file = open(W_FILE_PATH, 'r')

X = []
Y = []
Z = []

x_first = None

x_temp = []
y_temp = []
z_temp = []

for line in _w_file:
    line = line.replace('\n', '')
    line_splitted = line.split(',')

    x = float(line_splitted[0])
    y = float(line_splitted[1])
    z = float(line_splitted[2])

    if None != x_first and x == x_first:
        X.append(x_temp)
        Y.append(y_temp)
        Z.append(z_temp)

        x_temp = []
        y_temp = []
        z_temp = []

    if None == x_first:
        x_first = x

    x_temp.append(x)
    y_temp.append(y)
    z_temp.append(z)

_w_file.close()

figure = plt.figure()
axes = figure.gca(projection='3d')

surface = axes.plot_surface(np.asarray(X), np.asarray(Y), np.asarray(Z), rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
figure.colorbar(surface, shrink=0.5, aspect=5)

plt.savefig('test_4.png', dpi=1000)
#plt.show()
