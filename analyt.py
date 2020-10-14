import math
from utilities import plot_x_y
import matplotlib.pyplot as plt

# conc analyt
L = 10
conc_analyt = []
grid_block_n_analyt = 20
dX = L / grid_block_n_analyt

grid_centers_analyt = []
for i in range(grid_block_n_analyt):
    grid_centers_analyt.append(0 + i * dX + dX / 2)

conc_ini = 10.0
conc_out = 20.0

D = 3.E-1
t = 10.0
dt = 1.0

time = []

for i in range(grid_block_n_analyt):
    conc_it = conc_out + (conc_ini - conc_out) * math.erf(
        i * dX / 2 / math.sqrt(D * t))
    conc_analyt.append(conc_it)
    time.append(i * dX)

fig, ax = plt.subplots()

plot_x_y(ax, time, conc_analyt, 'time, sec', 'concentration, kg/m3', '-')
