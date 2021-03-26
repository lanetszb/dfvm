import math
from utilities import plot_x_y
import matplotlib.pyplot as plt

# conc analyt
L = 0.001
conc_analyt = []
grid_block_n_analyt = 100
dX = L / grid_block_n_analyt

grid_centers_analyt = []
for i in range(grid_block_n_analyt):
    grid_centers_analyt.append(0 + i * dX + dX / 2)

conc_ini = 1.0
conc_out = 10.0

D = 1.5E-7
t = 0.01

time = []

for i in range(grid_block_n_analyt):
    conc_it = conc_out + (conc_ini - conc_out) * math.erf(
        i * dX / 2 / math.sqrt(D * t))
    conc_analyt.append(conc_it)
    time.append(i * dX)

fig, ax = plt.subplots()

plot_x_y(ax, time, conc_analyt, 'time, sec', 'concentration, kg/m3', '-')
plt.show()
