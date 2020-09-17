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

#include "Equation.h"
#include "eigenSetGet.h"
#include <algorithm>
#include <iterator>
#include <numeric>

Equation::Equation(std::shared_ptr<Props> props,
                   std::shared_ptr<Sgrid> sgrid,
                   std::shared_ptr<Local> local,
                   std::shared_ptr<Convective> convective) :

        _props(props),
        _sgrid(sgrid),
        _local(local),
        _convective(convective),
        dim(_sgrid->_cellsN),
        iCurr(0), iPrev(1),
        _concsIni(new double[dim], dim),
        matrix(dim, dim),
        freeVector(new double[dim], dim) {

    std::vector<Triplet> triplets;
    triplets.reserve(3 * dim - 4);

    for (int i = 0; i < dim; i++)
        triplets.emplace_back(i, i);

    auto nonBoundCells = _sgrid->_typesCells.at("nonbound");
    for (int i = 0; i < nonBoundCells.size(); i++)
        for (int j = 0;
             j < _sgrid->_neighborsFaces.at(nonBoundCells[i]).size(); j++)
            for (int k = 0; k < _sgrid->_neighborsCells.at(
                    _sgrid->_neighborsFaces.at(
                            nonBoundCells[i])[j]).size(); k++) {

                auto neigborFace = _sgrid->_neighborsFaces.at(
                        nonBoundCells[i])[j];
                auto neigborCell = _sgrid->_neighborsCells.at(neigborFace)[k];

                if (neigborCell != nonBoundCells[i])
                    triplets.emplace_back(nonBoundCells[i], neigborCell);
            }

    matrix.setFromTriplets(triplets.begin(), triplets.end());
}

void Equation::fillMatrix(const double &alpha) {

    for (int i = 0; i < _sgrid->_cellsN; i++)
        matrix.coeffRef(i, i) = alpha;

    auto nonBoundCells = _sgrid->_typesCells.at("nonbound");
    for (int i = 0; i < nonBoundCells.size(); i++) {

        MatrixIterator matrixIt(matrix, nonBoundCells[i]);
        auto neighborFaces = _sgrid->_neighborsFaces.at(nonBoundCells[i]);
        auto facesIt = 0;

        while (facesIt < neighborFaces.size()) {
            if (matrixIt.col() != matrixIt.row()) {
                auto beta = _convective->_betas[neighborFaces[facesIt]];
                matrixIt.valueRef() = -1 * beta;
                matrix.coeffRef(nonBoundCells[i], nonBoundCells[i]) += beta;
                ++matrixIt;
                ++facesIt;
            } else {
                ++matrixIt;
            }
        }
    }
}

void Equation::setDirichletToBound(const double &concBound, const double &alpha,
                                   Eigen::Map<Eigen::VectorXi> boundCells) {

    for (int i = 0; i < boundCells.size(); i++)
        freeVector[boundCells[i]] = alpha * concBound;
}


void Equation::calculateFreeVector(const double &alpha) {

    auto concBack = 40.0;
    auto concBottom = 10.0;
    auto concFront = 20.0;
    auto concLeft = 30.0;
    auto concRight = 15.0;
    auto concTop = 27.0;

    auto cellsBack = _sgrid->_typesCells.at("back");
    auto cellsBottom = _sgrid->_typesCells.at("bottom");
    auto cellsFront = _sgrid->_typesCells.at("front");
    auto cellsLeft = _sgrid->_typesCells.at("left");
    auto cellsRight = _sgrid->_typesCells.at("right");
    auto cellsTop = _sgrid->_typesCells.at("top");

    setDirichletToBound(concBack, alpha, cellsBack);
    setDirichletToBound(concBottom, alpha, cellsBottom);
    setDirichletToBound(concFront, alpha, cellsFront);
    setDirichletToBound(concLeft, alpha, cellsLeft);
    setDirichletToBound(concRight, alpha, cellsRight);
    setDirichletToBound(concTop, alpha, cellsTop);

    auto nonBoundCells = _sgrid->_typesCells.at("nonbound");
    for (int i = 0; i < nonBoundCells.size(); i++)
        freeVector[nonBoundCells[i]] = alpha * _concs[iPrev][nonBoundCells[i]];
}

void Equation::calcConcsImplicit() {

    BiCGSTAB biCGSTAB;

    biCGSTAB.compute(matrix);

    _concs[iCurr] = biCGSTAB.solveWithGuess(freeVector, _concs[iPrev]);

}

void Equation::calcConcsExplicit() {

    for (int i = 0; i < _concsIni.size(); i++)
        _concsIni[i] = 20.0;
}

void Equation::cfdProcedure() {

    // ToDO: ask Alesk how to better do it
    _concs.emplace_back(_sgrid->_cellsArrays.at("concs_array1"));
    _concs.emplace_back(_sgrid->_cellsArrays.at("concs_array2"));

    _local->calcTimeSteps();
    _local->calcAlphas();

    for (int i = 0; i < _local->_alphas.size(); i++) {
        std::swap(iCurr, iPrev);
        _convective->calcBetas(_concs[iPrev]);
        fillMatrix(_local->_alphas[i]);
        calculateFreeVector(_local->_alphas[i]);
        calcConcsImplicit();
        Eigen::Map<Eigen::VectorXd> concCurr(new double[dim], dim);
        concCurr = _concs[iCurr];
        _concsTime.push_back(concCurr);
    }
}

std::vector<Eigen::Ref<Eigen::VectorXd>> Equation::getConcs() {
    return Eigen::vectorGetter<Eigen::VectorXd>(_concs);
}

void Equation::setConcs(std::vector<Eigen::Ref<Eigen::VectorXd>> &concs) {
    Eigen::vectorSetter<Eigen::VectorXd>(concs, _concs);
}

std::vector<Eigen::Ref<Eigen::VectorXd>> Equation::getConcsTime() {
    return Eigen::vectorGetter<Eigen::VectorXd>(_concsTime);
}

void
Equation::setConcsTime(std::vector<Eigen::Ref<Eigen::VectorXd>> &concsTime) {
    Eigen::vectorSetter<Eigen::VectorXd>(concsTime, _concsTime);
}

Eigen::Ref<Eigen::VectorXd> Equation::getConcsIni() {
    return _concsIni;
}

void Equation::setConcsIni(Eigen::Ref<Eigen::VectorXd> concsIni) {
    if (_concsIni.data() != concsIni.data())
        delete _concsIni.data();
    new(&_concsIni) Eigen::Map<Eigen::VectorXd>(concsIni.data(),
                                                concsIni.size());
}