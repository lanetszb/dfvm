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

import json

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))
print(current_path)

from diffusion_bind import Props
from diffusion_bind import Local
from diffusion_bind import Convective
from diffusion_bind import Equation
from sgrid import Sgrid

points_dims = np.array([5, 5, 5], dtype=np.int32)
points_origin = np.array([0., 0., 0.], dtype=np.float)
spacing = np.array([1., 1., 1.], dtype=np.float)

sgrid = Sgrid(points_dims, points_origin, spacing)

points_array = np.random.rand(sgrid.points_N)
points_arrays = {"points_array": points_array}
sgrid.points_arrays = points_arrays

cells_array = np.random.rand(sgrid.cells_N)
cells_arrays = {"cells_array ": cells_array}
sgrid.cells_arrays = cells_arrays

#
time_period = float(10)  # sec
time_step = 0.9  # sec
d_coeff_a = 0.0  # m2/sec
d_coeff_b = 1.E-10  # m2/sec
#
params = {'time_period': time_period, 'time_step': time_step,
          'd_coeff_a': d_coeff_a, 'd_coeff_b': d_coeff_b}
#
props = Props(params)
local = Local(props, sgrid)
convective = Convective(props, sgrid)

equation = Equation(props, sgrid, local, convective)


# local.calc_time_steps()
# local.calc_alphas()
#
# convective = Convective(props, sgrid)
# concs = np.ones(sgrid.cells_N)
# convective.calc_betas(concs)
# convective.weigh_D("mean_Average")
