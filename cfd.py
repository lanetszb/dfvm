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
import copy
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from dfvm import Props, Boundary, Local, Convective, Equation
from dfvm import calc_a_func, calc_b_func, calc_poro
from dfvm import plot_x_y
from sgrid import Sgrid, save_files_collection_to_file

# model geometry
points_dims = [5, 5, 2]
points_origin = [0., 0., 0.]
spacing = [1., 1., 1.]

sgrid = Sgrid(points_dims, points_origin, spacing)
points_array = np.random.rand(sgrid.points_N)

points_arrays = {"points_array": points_array}
active_cells = np.arange(sgrid.cells_N, dtype=np.uint64)
# initial concentration
conc_ini = float(0.0)
concs_array1 = np.tile(conc_ini, sgrid.cells_N)
concs_array2 = np.tile(conc_ini, sgrid.cells_N)
concs_arrays = {"concs_array1": concs_array1,
                "concs_array2": concs_array2}

sgrid.cells_arrays = concs_arrays
sgrid.set_cells_type('active', active_cells)
sgrid.process_type_by_cells_type('active')

# computation time
time_period = float(1000)  # sec
# numerical time step
time_step = float(1)  # sec

# diffusivity coeffs (specify only b coeff to make free diffusion constant)
d_coeff_a = float(0)  # m2/sec
d_coeff_b = float(15.E-3)  # m2/sec
# porosity of rock
poro_ini = float(1)
params = {'time_period': time_period, 'time_step': time_step,
          'd_coeff_a': d_coeff_a, 'd_coeff_b': d_coeff_b,
          'poro': poro_ini}

key_dirichlet_one = 'left'
key_dirichlet_two = 'right'

props = Props(params)
boundary = Boundary(props, sgrid)
boundary_faces_one = copy.deepcopy(sgrid.types_faces[key_dirichlet_one])
boundary_faces_two = copy.deepcopy(sgrid.types_faces[key_dirichlet_two])
boundary_face_one = sgrid.types_faces[key_dirichlet_one][0]
boundary_faces_one_axis = sgrid.faces_axes[boundary_face_one]
boundary_face_two = sgrid.types_faces[key_dirichlet_two][0]
boundary_faces_two_axis = sgrid.faces_axes[boundary_face_two]
boundary.shift_boundary_faces(boundary_faces_one, boundary_faces_one_axis)
boundary.shift_boundary_faces(boundary_faces_two, boundary_faces_two_axis)

local = Local(props, sgrid)
convective = Convective(props, sgrid)
equation = Equation(props, sgrid, local, convective)

# dirichlet cells (options: left, right, top, bottom, front, back)
equation.bound_groups_dirich = [key_dirichlet_one, key_dirichlet_two]
# concentration on dirichlet cells
conc_left = float(20)
conc_right = float(0)
equation.concs_bound_dirich = {key_dirichlet_one: conc_left,
                               key_dirichlet_two: conc_right}

# equation.cfd_procedure()

# CFD procedure
concs = [concs_array1, concs_array2]
equation.concs = concs

local.calc_time_steps()
time_steps = local.time_steps
concs_time = []
conc_curr = copy.deepcopy(equation.concs[equation.i_curr])
concs_time.append(conc_curr)
flow_rate_one_time = []
flow_rate_two_time = []
for time_step in time_steps:
    # modifing porosity
    equation.cfd_procedure_one_step(time_step)
    conc_curr = copy.deepcopy(equation.concs[equation.i_curr])
    concs_time.append(conc_curr)
    flow_rate_boundary_one = equation.calc_faces_flow_rate(boundary_faces_one)
    flow_rate_boundary_two = equation.calc_faces_flow_rate(boundary_faces_two)
    flow_rate_one_time.append(flow_rate_boundary_one)
    flow_rate_two_time.append(flow_rate_boundary_two)
    # new Dirichlet boundaries can be input here
    equation.concs_bound_dirich = {key_dirichlet_one: conc_left, key_dirichlet_two: conc_right}
equation.concs_time = concs_time
#

# visualising 'a' and 'b' coefficients and porosity
# set concentration range for visualisation
a_list = []
b_list = []
poro_list = []
conc_list = []
for i in range(int(conc_right), int(conc_left)):
    conc_list.append(i)
    a_list.append(calc_a_func(i, poro_ini))
    b_list.append(calc_b_func(i, d_coeff_b, poro_ini))
    poro_list.append(calc_poro(i, poro_ini))

# plotting the dependence of 'a' and 'b' coefficients and porosity on free concentration
fig, axs = plt.subplots(3, sharex=True)
plot_x_y(axs[0], conc_list, a_list, 'concentration', 'coeff a', '-',
         color='green')
plot_x_y(axs[1], conc_list, b_list, 'concentration', 'coeff b', '-',
         color='blue')
plot_x_y(axs[2], conc_list, poro_list, 'concentration', 'poro', '-',
         color='red')
axs[0].legend('a', loc="best")
axs[1].legend('b', loc="best")
axs[2].legend(['poro'], loc="best")

# plotting inlet and outlet flow rates vs time
time = np.cumsum(np.array(time_steps))
fig1, ax1 = plt.subplots()
plot_x_y(ax1, time, flow_rate_one_time, 'time', 'G, kg/sec', '-',
         color='green')
plot_x_y(ax1, time, flow_rate_two_time, 'time', 'G, kg/sec', '-',
         color='blue')
ax1.legend(['$Q_{in}$', '$Q_{out}$'], loc="best")


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
    sgrid.save_cells('inOut/' + files_names[i])

save_files_collection_to_file(file_name, files_names, files_descriptions)
