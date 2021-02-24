from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

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

axes.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True)
plt.show()
