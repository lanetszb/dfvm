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
        freeVector(new double[dim], dim) {}


void Equation::procesNewmanFaces(const double &flowNewman,
                                 Eigen::Map<Eigen::VectorXui64> faces) {

    for (int i = 0; i < faces.size(); i++) {
        auto &face = faces[i];
        auto &cells = _sgrid->_neighborsCells.at(face);
        for (int j = 0; j < cells.size(); j++) {
            auto &cell = cells[j];
            _coeffsMatrix[face][cell] = 0;
            _coeffsFreeVec[face][cell] = flowNewman;
        }
    }

}

// ToDo Ref instead of Map
void Equation::processNonBoundFaces(Eigen::Map<Eigen::VectorXui64> faces) {

    for (int i = 0; i < faces.size(); i++) {
        auto &face = faces[i];
        auto &cells = _sgrid->_neighborsCells.at(face);
        for (int j = 0; j < cells.size(); j++) {
            auto &cell = cells[j];
            auto &normal = _sgrid->_normalsNeighborsCells.at(face)[j];
            _coeffsMatrix[face][cell] = normal * _convective->_betas[face];
            _coeffsFreeVec[face][cell] = 0;
        }
    }

}

void Equation::fillMatrix(const double &alpha) {

    matrix.setZero();

    std::vector<std::string> dirichCellsBounds = {"left", "right",
                                                  "top", "bottom",
                                                  "front", "back"};
    std::vector<int> dirichCells;
    for (auto &bound : dirichCellsBounds) {
        auto cells = _sgrid->_typesCells.at(bound);
        auto position = dirichCells.end();
        dirichCells.insert(position, cells.data(),
                           cells.data() + cells.size());
    }

    // std::vector<std::string> dirichCellsBounds = {"left", "right"};

    // std::vector<int> dirichCells;
    // for (auto &bound : dirichCellsBounds) {
    //     auto cells = _sgrid->_typesCells.at(bound);
    //     auto position = dirichCells.end();
    //     dirichCells.insert(position, cells.data(),
    //                        cells.data() + cells.size());
    // }

    for (auto &cell: dirichCells)
        matrix.coeffRef(cell, cell) = alpha;

    auto nonBoundCells = _sgrid->_typesCells.at("nonbound");

    for (int i = 0; i < nonBoundCells.size(); i++)
        matrix.coeffRef(nonBoundCells[i], nonBoundCells[i]) = alpha;

    for (int i = 0; i < nonBoundCells.size(); i++) {
        auto faces = _sgrid->_neighborsFaces.at(nonBoundCells[i]);
        auto normalsFaces = _sgrid->_normalsNeighborsFaces.at(nonBoundCells[i]);
        for (int j = 0; j < faces.size(); j++)
            for (const auto&[cell, cellCoeff] : _coeffsMatrix[faces[j]])
                matrix.coeffRef(nonBoundCells[i], cell) +=
                        normalsFaces[j] * cellCoeff;
    }

    // std::vector<std::string> noFlowCellsBounds = {"top", "bottom",
    //                                               "front", "back"};

    // std::vector<int> noFlowCells;
    // for (auto &bound : noFlowCellsBounds) {
    //     auto cells = _sgrid->_typesCells.at(bound);
    //     auto position = noFlowCells.end();
    //     noFlowCells.insert(position, cells.data(),
    //                        cells.data() + cells.size());
    // }
    // sort(dirichCells.begin(), dirichCells.end());
    // sort(noFlowCells.begin(), noFlowCells.end());
    // noFlowCells.erase(unique(noFlowCells.begin(), noFlowCells.end()),
    //                   noFlowCells.end());

    // std::vector<int> diff;

    // std::set_difference(std::begin(noFlowCells), std::end(noFlowCells),
    //                     std::begin(dirichCells), std::end(dirichCells),
    //                     std::back_inserter(diff));

    // for (const auto &i: dirichCells)
    //     std::cout << i << ' ';
    // std::cout << std::endl;

    // for (const auto &i: noFlowCells)
    //     std::cout << i << ' ';
    // std::cout << std::endl;

    // for (const auto &i: diff)
    //     std::cout << i << ' ';
    // std::cout << std::endl;

    // for (int i = 0; i < diff.size(); i++)
    //     matrix.coeffRef(diff[i], diff[i]) = alpha;

    // for (int i = 0; i < diff.size(); i++) {
    //     auto faces = _sgrid->_neighborsFaces.at(diff[i]);
    //     auto normalsFaces = _sgrid->_normalsNeighborsFaces.at(diff[i]);
    //     for (int j = 0; j < faces.size(); j++)
    //         for (const auto&[cell, cellCoeff] : _coeffsMatrix[faces[j]])
    //             matrix.coeffRef(diff[i], cell) +=
    //                     normalsFaces[j] * cellCoeff;
    // }
}

void Equation::procesDirichCells(const double &alpha,
                                 std::vector<std::string> &boundGroups,
                                 std::map<std::string, double> &concsBound) {


    for (int i = 0; i < boundGroups.size(); i++) {
        auto &bound = boundGroups[i];
        auto &conc = concsBound[bound];
        auto cells = _sgrid->_typesCells.at(bound);
        for (int j = 0; j < cells.size(); j++) {
            auto cell = cells[j];
            freeVector[cell] = conc * alpha;
        }
    }

}

void Equation::calculateFreeVector(const double &alpha) {

    procesDirichCells(alpha, _boundGroupsDirich, _concsBoundDirich);

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

        // procesNoFlowFaces(_sgrid->_typesFaces.at("top"));
        // procesNoFlowFaces(_sgrid->_typesFaces.at("bottom"));
        // procesNoFlowFaces(_sgrid->_typesFaces.at("front"));
        // procesNoFlowFaces(_sgrid->_typesFaces.at("back"));
        processNonBoundFaces(_sgrid->_typesFaces.at("nonbound"));

        fillMatrix(_local->_alphas[i]);
        calculateFreeVector(_local->_alphas[i]);
        calcConcsImplicit();
        Eigen::Map<Eigen::VectorXd> concCurr(new double[dim], dim);
        concCurr = _concs[iCurr];
        _concsTime.push_back(concCurr);

        // std::cout << matrix << std::endl;
    }

    // for (const auto&[face, faceN] : _coeffsMatrix)
    //     for (const auto&[cell, cellCoeff] : faceN)
    //         std::cout << "coeffs" << "[" << face << "][" << cell << "] = "
    //                   << cellCoeff
    //                   << std::endl;
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