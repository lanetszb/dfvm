/* MIT License
 *
 * Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "math/Props.h"
#include "math/Local.h"
#include "math/Convective.h"
// #include "Equation.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(diffusion_bind, m) {

    py::class_<Props, std::shared_ptr<Props>>(m, "Props")
            .def(py::init<const std::map<std::string, std::variant<int, double>> &>(),
                 "params"_a)

            .def_readwrite("params", &Props::_params)
            .def("calc_D", &Props::calcD, "conc"_a)
            .def("print_params", &Props::printParams);

    py::class_<Local, std::shared_ptr<Local>>(m, "Local")
            .def(py::init<std::shared_ptr<Props>, std::shared_ptr<Sgrid>>(),
                 "props"_a,
                 "sgrid"_a)

            .def("calc_time_steps", &Local::calcTimeSteps)
            .def("calc_alphas", &Local::calcAlphas)
            .def_property("time_steps",
                          &Local::getTimeSteps, &Local::setTimeSteps)
            .def_property("alphas",
                          &Local::getAlphas, &Local::setAlphas);


    py::class_<Convective, std::shared_ptr<Convective>>(m, "Convective")
            .def(py::init<std::shared_ptr<Props>, std::shared_ptr<Sgrid>>(),
                 "props"_a,
                 "sgrid"_a)

            .def("weigh_D", &Convective::weighD, "method"_a)
            .def("calc_betas", &Convective::calcBetas)
            .def_property("betas",
                          &Convective::getBetas, &Convective::setBetas);

}


//            .def_readwrite("params", &Props::_params)
//            .def("calc_d", &Props::calcD, "conc"_a)
//            .def("print_params", &Props::printParams);


//  py::class_<Equation>(m, "Equation")
//      .def(py::init<const std::map<std::string, std::variant<int, double>> &>(),
//           "params"_a)
//
//      .def("calculate", &Equation::calculate)
//
//      .def_readwrite("concs", &Equation::concs)
//      .def_readwrite("velocities", &Equation::velocities)
//      .def_readwrite("time_steps", &Equation::_timeSteps)
//      .def_readwrite("times", &Equation::times);
//
//}