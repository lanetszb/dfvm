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

#include "Local.h"
#include <iostream>
#include <algorithm>
#include <iterator>

Local::Local(std::shared_ptr<Props> props, std::shared_ptr<Sgrid> sgrid) :
        _props(props),
        _sgrid(sgrid),
        _alphas(new double, 1),
        _timeSteps(new double, 1) {
}

void Local::calculateTimeSteps() {

    auto time = std::get<double>(_props->_params["time_period"]);
    auto timeStep = std::get<double>(_props->_params["time_step"]);
    double division = time / timeStep;
    double fullStepsN;
    auto lastStep = std::modf(division, &fullStepsN);
    auto localTimeSteps = std::vector<double>(fullStepsN, timeStep);
    if (lastStep > 0)
        localTimeSteps.push_back(lastStep * timeStep);

    auto array = new double[localTimeSteps.size()];
    std::copy(std::begin(localTimeSteps), std::end(localTimeSteps), array);

    delete _timeSteps.data();
    new(&_timeSteps) Eigen::Map<Eigen::VectorXd>(array, localTimeSteps.size());
    std::cout << _timeSteps << std::endl;
}

void Local::calculateAlphas() {

    delete _alphas.data();
    auto array = new double[_timeSteps.size()];

    new(&_alphas) Eigen::Map<Eigen::VectorXd>(array, _timeSteps.size());

    _alphas = _sgrid->_cellV * _timeSteps.cwiseInverse();
}

Eigen::Ref<Eigen::VectorXd> Local::getTimeSteps() {
    return _timeSteps;
}

void Local::setTimeSteps(Eigen::Ref<Eigen::VectorXd> timeSteps) {
    new(&_timeSteps) Eigen::Map<Eigen::VectorXd>(timeSteps.data(),
                                                 timeSteps.size());
}

Eigen::Ref<Eigen::VectorXd> Local::getAlphas() {
    return _alphas;
}

void Local::setAlphas(Eigen::Ref<Eigen::VectorXd> alphas) {
    new(&_alphas) Eigen::Map<Eigen::VectorXd>(alphas.data(),
                                              alphas.size());
}