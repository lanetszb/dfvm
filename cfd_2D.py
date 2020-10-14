# MIT License
#
# Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import sys
import os
import numpy as np
import matplotlib.pyplot as plt

from utilities import plot_x_y

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))
print(current_path)

import diffusion_bind
from diffusion_bind import Props, Local, Convective, Equation
from sgrid import Sgrid, save_files_collection_to_file

# model geometry
points_dims = [201, 201, 2]
points_origin = [0., 0., 0.]
spacing = [1., 1., 1.]

sgrid = Sgrid(points_dims, points_origin, spacing)

# initial concentration
conc_ini = 10.
concs_array1 = np.tile(conc_ini, sgrid.cells_N)
concs_array2 = np.tile(conc_ini, sgrid.cells_N)
concs_arrays = {"concs_array1": concs_array1,
                "concs_array2": concs_array2}
sgrid.cells_arrays = concs_arrays
# computation time
time_period = 50.  # sec
# numerical time step
time_step = 5.  # sec

# diffusivity coeffs (specify only b coeff to make free diffusion constant)
d_coeff_a = float(0.0)  # m2/sec
d_coeff_b = float(3.E-1)  # m2/sec
# porosity of rock
poro = float(1)
params = {'time_period': time_period, 'time_step': time_step,
          'd_coeff_a': d_coeff_a, 'd_coeff_b': d_coeff_b,
          'poro': poro}

props = Props(params)
local = Local(props, sgrid)
convective = Convective(props, sgrid)
equation = Equation(props, sgrid, local, convective)

# dirichlet cells (options: left, right, top, bottom, front, back)
equation.bound_groups_dirich = ['left', 'right']
# concentration on dirichlet cells
conc_left = float(10)
conc_right = float(20)
equation.concs_bound_dirich = {'left': conc_left, 'right': conc_right}
equation.cfd_procedure()

# saving results to paraview

os.system('rm -r inOut/*.vtu')
os.system('rm -r inOut/*.pvd')
concs_dict = dict()
file_name = 'inOut/collection.pvd'
files_names = list()
files_descriptions = list()
for i in range(len(local.time_steps)):
    sgrid.cells_arrays = {'conc_i': equation.concs_time[i]}
    files_names.append(str(i) + '.vtu')
    files_descriptions.append(str(i))
    sgrid.save('inOut/' + files_names[i])

save_files_collection_to_file(file_name, files_names, files_descriptions)
