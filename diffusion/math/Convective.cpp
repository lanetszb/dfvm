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
#include <algorithm>
#include <iterator>

Convective::Convective(std::shared_ptr<Props> props,
                       std::shared_ptr<Sgrid> sgrid) :
        _props(props),
        _sgrid(sgrid),
        _betas(_sgrid->_facesN) {}

double Convective::weighD(const std::string &method, const double &conc0,
                          const double &conc1) {

    if (method == "meanAverage") {
        auto conc = (conc0 + conc1) / 2;
        return _props->calcD(conc);

    } else if (method == "meanHarmonic") {
        auto conc = 2. * conc0 * conc1 / (conc0 + conc1);
        return _props->calcD(conc);

    } else if (method == "upWind") {
        auto conc = std::max(conc0, conc1);
        return _props->calcD(conc);

    } else exit(0);
}

void Convective::calcBetas(Eigen::Ref<Eigen::VectorXd> concs, double &time) {

    clock_t tStart = clock();

    auto &neighborsCells = _sgrid->_neighborsCells;

    for (int i = 0; i < _betas.size(); i++) {
        auto &neighborsCell = neighborsCells[i];
        auto &conc0 = concs(neighborsCell[0]);
        auto &conc1 = concs(neighborsCell[1]);

        double diffusivity;
        if (neighborsCells[i].size() == 2) {
            diffusivity = weighD("meanAverage", conc0, conc1);
        }
        else {
            diffusivity = weighD("meanAverage", conc0, conc0);
        }

        auto &axis = _sgrid->_facesAxes[i];
        _betas[i] =
                diffusivity * _sgrid->_facesSs[axis] / _sgrid->_spacing[axis];

    }
    time += (double) (clock() - tStart);
}