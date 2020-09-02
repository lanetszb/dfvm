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
        _betas(new double[_sgrid->_facesN], _sgrid->_facesN) {}

double Convective::weighD(const std::string &method) {

    if (method == "meanAverage")
        return 1.0;
    else if (method == "meanHarmonic")
        return 2.0;
    else if (method == "upWind")
        return 3.0;
    else exit(0);
}

void Convective::calcBetas(Eigen::Ref<Eigen::VectorXd> concs) {

    _betas = concs;
}

Eigen::Ref<Eigen::VectorXd> Convective::getBetas() {
    return _betas;
}

void Convective::setBetas(Eigen::Ref<Eigen::VectorXd> betas) {
    if (_betas.data() != betas.data())
        delete _betas.data();
    new(&_betas) Eigen::Map<Eigen::VectorXd>(betas.data(),
                                             betas.size());
}