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

#include "Convective.h"
#include "funcs.h"
#include <algorithm>
#include <iterator>

Convective::Convective(std::shared_ptr<Props> props,
                       std::shared_ptr<Sgrid> sgrid) :
        _props(props),
        _sgrid(sgrid),
        _betas(_sgrid->_facesN) {}
// ToDo: weighing (method, valueA, valueB)
double Convective::weighBcoeff(const std::string &method, const double &conc0,
                               const double &conc1) {

    auto poro = std::get<double>(_props->_params["poro"]);
    auto diffusivity0 = _props->calcD(conc0);
    auto bCoeff0 = calcBFunc(conc0, diffusivity0, poro);

    auto diffusivity1 = _props->calcD(conc1);
    auto bCoeff1 = calcBFunc(conc1, diffusivity1, poro);

    if (method == "meanAverage")
        return (bCoeff0 + bCoeff1) / 2;
    else if (method == "meanHarmonic")
        return 2. * bCoeff0 * bCoeff1 / (bCoeff0 + bCoeff1);
    else if (method == "upWind")
        return std::max(bCoeff0, bCoeff1);
    else exit(0);
}

// ToDo: massive of diffusions which are going to be different for matrix and fractures
void Convective::calcBetas(Eigen::Ref<Eigen::VectorXd> concs, double &time) {

    clock_t tStart = clock();

    auto &neighborsCells = _sgrid->_neighborsCells;
    auto boundFaces = _sgrid->_typesFaces.at("active_bound");
    auto nonBoundFaces = _sgrid->_typesFaces.at("active_nonbound");

    for (int i = 0; i < boundFaces.size(); i++) {
        auto boundFace = boundFaces[i];
        auto &faceNeighborsCell = neighborsCells[boundFace];
        auto &conc0 = concs(faceNeighborsCell[0]);

        auto bCoeff = weighBcoeff("meanAverage", conc0, conc0);
        auto &axis = _sgrid->_facesAxes[boundFace];

        _betas[boundFace] =
                bCoeff * _sgrid->_facesSs[axis] / _sgrid->_spacing[axis];
    }

    for (int i = 0; i < nonBoundFaces.size(); i++) {
        auto nonBoundFace = nonBoundFaces[i];
        auto &faceNeighborsCell = neighborsCells[nonBoundFace];
        auto &conc0 = concs(faceNeighborsCell[0]);
        auto &conc1 = concs(faceNeighborsCell[1]);

        auto bCoeff = weighBcoeff("meanAverage", conc0, conc1);
        auto &axis = _sgrid->_facesAxes[nonBoundFace];

        _betas[nonBoundFace] = bCoeff * _sgrid->_facesSs[axis]
                               / _sgrid->_spacing[axis];
    }

    time += (double) (clock() - tStart);
}