from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import matplotlib.pyplot as plt
import numpy as np

W_FILE_PATH = '/home/glen/Projects/ice_deflection/data/w_values.txt'

_w_file = open(W_FILE_PATH, 'r')

x = []
y = []
z = []

for line in _w_file:
    line = line.replace('\n', '')
    line_splitted = line.split(',')

    x.append(float(line_splitted[0]))
    y.append(float(line_splitted[1]))
    z.append(float(line_splitted[2]))

_w_file.close()

figure = plt.figure()
axes = figure.gca(projection='3d')

surface = axes.plot_trisurf(x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
figure.colorbar(surface, shrink=0.5, aspect=5)

plt.savefig('data/test.png', dpi=1000)
plt.show()
