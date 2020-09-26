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
#include <time.h>
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

    std::vector<std::string> boundGroupsDirich{"left", "right"};
    auto nonBoundCells = findNonDirichCells(boundGroupsDirich);

    for (int i = 0; i < nonBoundCells.size(); i++) {
        auto faces = _sgrid->_neighborsFaces.at(nonBoundCells[i]);
        for (int j = 0; j < faces.size(); j++) {
            auto &face = faces[j];
            auto cells = _sgrid->_neighborsCells.at(face);
            for (int k = 0; k < cells.size(); k++) {
                auto &cell = cells[k];
                triplets.emplace_back(nonBoundCells[i], cell);
            }

        }
    }

    matrix.setFromTriplets(triplets.begin(), triplets.end());
}


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
void Equation::processNonBoundFaces(Eigen::Ref<Eigen::VectorXui64> faces) {

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

std::vector<int>
Equation::groupVecsByKeys(std::vector<std::string> &groups) {

    std::vector<int> groupedCells;
    for (auto &bound : groups) {
        auto cells = _sgrid->_typesCells.at(bound);
        auto position = groupedCells.end();
        groupedCells.insert(position, cells.data(),
                            cells.data() + cells.size());
    }

    return groupedCells;
}


std::vector<int> Equation::findNonDirichCells
        (std::vector<std::string> &boundGroupsDirich) {

    std::vector<int> dirichCells = groupVecsByKeys(boundGroupsDirich);

    std::vector<std::string> groupsNonDirich;
    auto allCells = _sgrid->_typesCells;
    for (auto &type : allCells)
        for (auto &typeDirich : boundGroupsDirich) {
            if (typeDirich != type.first)
                groupsNonDirich.push_back(type.first);
        }

    std::vector<int> nonDirichCellsDump = groupVecsByKeys(groupsNonDirich);

    sort(dirichCells.begin(), dirichCells.end());
    sort(nonDirichCellsDump.begin(), nonDirichCellsDump.end());
    nonDirichCellsDump.erase(
            unique(nonDirichCellsDump.begin(), nonDirichCellsDump.end()),
            nonDirichCellsDump.end());

    std::vector<int> nonDirichCells;

    std::set_difference(std::begin(nonDirichCellsDump),
                        std::end(nonDirichCellsDump),
                        std::begin(dirichCells), std::end(dirichCells),
                        std::back_inserter(nonDirichCells));

    return nonDirichCells;
}

void Equation::fillMatrix(const double &alpha,
                          std::map<std::string, double> &times) {
    clock_t tStart = clock();

    for (int i = 1; i < dim; ++i) {
        for (MatrixIterator it(matrix, i); it; ++it)
            it.valueRef() = 0;
    }

    auto &cellsN = _sgrid->_cellsN;
    for (uint64_t i = 0; i < cellsN; i++) {
        matrix.coeffRef(i, i) = alpha;
        freeVector[i] = alpha * _concs[iPrev][i];
    }

    auto nonBoundCells = findNonDirichCells(_boundGroupsDirich);
    for (int i = 0; i < nonBoundCells.size(); i++) {
        clock_t wholeLoopTimeStart = clock();
        auto faces = _sgrid->_neighborsFaces.at(nonBoundCells[i]);
        auto normalsFaces = _sgrid->_normalsNeighborsFaces.at(
                nonBoundCells[i]);
        for (int j = 0; j < faces.size(); j++) {
            clock_t innerLoopTimeStart = clock();
            for (const auto&[cell, cellCoeff] : _coeffsMatrix[faces[j]]) {
                clock_t innerInnerLoopTimeStart = clock();
                matrix.coeffRef(nonBoundCells[i], cell) +=
                        normalsFaces[j] * cellCoeff;
                times["loop1LoopLoop"] += (double) (clock() -
                                                    innerInnerLoopTimeStart);
            }
            times["loop1Loop"] += (double) (clock() - innerLoopTimeStart);
        }
        times["loop1"] += (double) (clock() - wholeLoopTimeStart);
    }

    times["whole"] += (double) (clock() - tStart);
    // std::cout << matrix << std::endl;
}

void Equation::procesDirichCells(const double &alpha,
                                 std::vector<std::string> &boundGroups,
                                 std::map<std::string, double> &concsBound,
                                 double &time) {

    clock_t tStart = clock();
    for (int i = 0; i < boundGroups.size(); i++) {
        auto &bound = boundGroups[i];
        auto &conc = concsBound[bound];
        auto cells = _sgrid->_typesCells.at(bound);
        for (int j = 0; j < cells.size(); j++) {
            auto cell = cells[j];
            freeVector[cell] = conc * alpha;
        }
    }
    time += (double) (clock() - tStart);
}

void Equation::calcConcsImplicit(double &time) {

    clock_t tStart = clock();
    BiCGSTAB biCGSTAB;

    biCGSTAB.compute(matrix);

    _concs[iCurr] = biCGSTAB.solveWithGuess(freeVector, _concs[iPrev]);
    time += (double) (clock() - tStart);

}

void Equation::calcConcsExplicit() {

    for (int i = 0; i < _concsIni.size(); i++)
        _concsIni[i] = 20.0;
}

void Equation::cfdProcedure() {

    _concs.emplace_back(_sgrid->_cellsArrays.at("concs_array1"));
    _concs.emplace_back(_sgrid->_cellsArrays.at("concs_array2"));

    _local->calcTimeSteps();
    _local->calcAlphas();

    double calcBetasTime = 0;
    std::map<std::string, double> fillMatrixTimes;
    fillMatrixTimes["whole"] = 0;
    fillMatrixTimes["loop1"] = 0;
    fillMatrixTimes["loop1Loop"] = 0;
    fillMatrixTimes["loop1LoopLoop"] = 0;
    double procesDirichCellsTime = 0;
    double calcConcsImplicitTime = 0;
    double totalTime = 0;

    auto alphas = _local->_alphas;
    for (int i = 0; i < alphas.size(); i++) {
        clock_t tStart = clock();
        auto &alpha = alphas[i];
        std::swap(iCurr, iPrev);
        _convective->calcBetas(_concs[iPrev], calcBetasTime);

        processNonBoundFaces(_sgrid->_typesFaces.at("nonbound"));
        fillMatrix(alpha, fillMatrixTimes);
        procesDirichCells(alpha, _boundGroupsDirich, _concsBoundDirich,
                          procesDirichCellsTime);

        calcConcsImplicit(calcConcsImplicitTime);
        Eigen::Map<Eigen::VectorXd> concCurr(new double[dim], dim);
        concCurr = _concs[iCurr];
        _concsTime.push_back(concCurr);

        totalTime += (double) (clock() - tStart);

    }
    printf("calcBetasTime: %.2fs\n", calcBetasTime / CLOCKS_PER_SEC);
    for (auto const&[type, time] : fillMatrixTimes) {
        std::cout << type;
        printf(": %.2fs\n", time / CLOCKS_PER_SEC);
    }


    printf("procesDirichCellsTime: %.2fs\n",
           procesDirichCellsTime / CLOCKS_PER_SEC);
    printf("calcConcsImplicitTime: %.2fs\n",
           calcConcsImplicitTime / CLOCKS_PER_SEC);
    printf("totalTime: %.2fs\n", totalTime / CLOCKS_PER_SEC);
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
Equation::setConcsTime(
        std::vector<Eigen::Ref<Eigen::VectorXd>> &concsTime) {
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