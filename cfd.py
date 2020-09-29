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
import math

import json

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))
print(current_path)

from diffusion_bind import Props, Local, Convective, Equation
from sgrid import Sgrid, save_files_collection_to_file

points_dims = [10, 10, 10]
points_origin = [0., 0., 0.]
spacing = [1., 1., 1.]

sgrid = Sgrid(points_dims, points_origin, spacing)
points_array = np.random.rand(sgrid.points_N)

points_arrays = {"points_array": points_array}
# cells_array = np.random.rand(sgrid.cells_N)
# cells_arrays = {"cells_array ": cells_array}
# sgrid.cells_arrays = cells_arrays

concs_array1 = np.tile(float(10.0), sgrid.cells_N)
concs_array2 = np.tile(float(10.0), sgrid.cells_N)
concs_arrays = {"concs_array1": concs_array1,
                "concs_array2": concs_array2}

sgrid.cells_arrays = concs_arrays
time_period = float(100.)  # sec
time_step = 5.0  # sec
d_coeff_a = 0.0  # m2/sec
d_coeff_b = 3.E-1  # m2/sec
#
params = {'time_period': time_period, 'time_step': time_step,
          'd_coeff_a': d_coeff_a, 'd_coeff_b': d_coeff_b}
#
props = Props(params)
local = Local(props, sgrid)
convective = Convective(props, sgrid)

equation = Equation(props, sgrid, local, convective)
equation.bound_groups_dirich = ['left', 'right']
equation.concs_bound_dirich = {'left': float(10), 'right': float(20)}
equation.cfd_procedure()
# equation.bound_groups_dirich = ['left', 'right',
#                                 'top', 'bottom',
#                                 'front', 'back']
#
# equation.concs_bound_dirich = {'left': float(10), 'right': float(20),
#                                'top': float(30), 'bottom': float(40),
#                                'front': float(50), 'back': float(60)}


os.system('rm -r inOut/*.vtu')
os.system('rm -r inOut/*.pvd')
concs_dict = dict()
file_name = 'inOut/collection.pvd'
files_names = list()
files_descriptions = list()
for i in range(len(local.alphas)):
    sgrid.cells_arrays = {'conc': equation.concs_time[i]}
    files_names.append(str(i) + '.vtu')
    files_descriptions.append(str(i))
    sgrid.save('inOut/' + files_names[i])

save_files_collection_to_file(file_name, files_names, files_descriptions)

# conc analyt
# L = 10
# conc_analyt = []
# grid_block_n_analyt = 20
# dX = L / grid_block_n_analyt
# # #
# grid_centers_analyt = []
# for i in range(grid_block_n_analyt):
#     grid_centers_analyt.append(0 + i * dX + dX / 2)
#
# conc_ini = 10.0
# conc_out = 20.0
# #
# D = d_coeff_b
# t = time_period
# dt = time_step
# # #
# for i in range(grid_block_n_analyt):
#     conc_it = conc_out + (conc_ini - conc_out) * math.erf(
#         i * dX / 2 / math.sqrt(D * t))
#     conc_analyt.append(conc_it)
# # #
# conc_analyt.reverse()
#

# local.calc_time_steps()
# local.calc_alphas()
#
# convective = Convective(props, sgrid)
# concs = np.ones(sgrid.cells_N)
# convective.calc_betas(concs)
# convective.weigh_D("mean_Average")
