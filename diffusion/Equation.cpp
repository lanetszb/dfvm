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
#include "math/funcs.h"

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

    for (auto &nonDirichCell: findNonDirichCells(_boundGroupsDirich))
        for (auto &face: _sgrid->_neighborsFaces[nonDirichCell])
            for (auto &cell: _sgrid->_neighborsCells[face])
                triplets.emplace_back(nonDirichCell, cell);

    matrix.setFromTriplets(triplets.begin(), triplets.end());
}


void Equation::processNewmanFaces(const double &flowNewman,
                                  Eigen::Map<Eigen::VectorXui64> faces) {

    for (int i = 0; i < faces.size(); i++) {
        auto &face = faces[i];
        for (auto &cell : _sgrid->_neighborsCells[face]) {
            _matrixFacesCells[face][cell] = 0;
            _freeFacesCells[face][cell] = flowNewman;
        }
    }

}

void Equation::processNonBoundFaces(Eigen::Ref<Eigen::VectorXui64> faces) {

    for (int i = 0; i < faces.size(); i++) {
        auto &face = faces[i];
        auto &cells = _sgrid->_neighborsCells[face];
        for (int j = 0; j < cells.size(); j++) {
            auto &cell = cells[j];
            auto &normal = _sgrid->_normalsNeighborsCells[face][j];
            _matrixFacesCells[face][cell] = normal * _convective->_betas[face];
            _freeFacesCells[face][cell] = 0;
        }
    }

}

std::vector<uint64_t> Equation::groupCellsByTypes
        (const std::vector<std::string> &groups) {

    std::vector<uint64_t> groupedCells;
    for (auto &bound : groups) {
        auto cells = _sgrid->_typesCells.at(bound);
        auto position = groupedCells.end();
        groupedCells.insert(position, cells.data(),
                            cells.data() + cells.size());
    }

    sort(groupedCells.begin(), groupedCells.end());
    groupedCells.erase(unique(groupedCells.begin(), groupedCells.end()),
                       groupedCells.end());

    return groupedCells;
}

std::vector<uint64_t> Equation::findNonDirichCells(
        std::vector<std::string> &boundGroupsDirich) {

    auto dirichCells = groupCellsByTypes(boundGroupsDirich);
    auto activeCells = groupCellsByTypes({"active"});

    std::vector<uint64_t> nonDirichCells;
    std::set_difference(std::begin(activeCells),
                        std::end(activeCells),
                        std::begin(dirichCells), std::end(dirichCells),
                        std::back_inserter(nonDirichCells));

    return nonDirichCells;
}

void Equation::fillMatrix() {

    for (int i = 0; i < dim; ++i)
        for (MatrixIterator it(matrix, i); it; ++it)
            it.valueRef() = 0;

    for (uint64_t i = 0; i < _sgrid->_cellsN; i++) {
        matrix.coeffRef(i, i) = _local->_alphas[i];
        freeVector[i] = _local->_alphas[i] * _concs[iPrev][i];
    }

    for (auto &nonDirichCell: findNonDirichCells(_boundGroupsDirich)) {

        auto faces = _sgrid->_neighborsFaces[nonDirichCell];
        auto normalsFaces = _sgrid->_normalsNeighborsFaces[nonDirichCell];
        for (int j = 0; j < faces.size(); j++) {
            auto face = faces[j];
            auto normalFace = normalsFaces[j];
            for (const auto&[cell, cellCoeff] : _matrixFacesCells[face])
                matrix.coeffRef(nonDirichCell, cell) += normalFace * cellCoeff;
        }
    }

}

void Equation::processDirichCells(std::vector<std::string> &boundGroups,
                                  std::map<std::string, double> &concsBound) {

    auto activeBoundCells = groupCellsByTypes({"active_bound"});

    for (auto &bound : boundGroups) {
        auto &conc = concsBound[bound];
        auto dirichCells = groupCellsByTypes({bound});

        std::vector<int> dirichCellsActive;

        set_intersection(activeBoundCells.begin(),
                         activeBoundCells.end(),
                         dirichCells.begin(),
                         dirichCells.end(),
                         std::back_inserter(dirichCellsActive));

        for (int j = 0; j < dirichCellsActive.size(); j++) {
            auto cell = dirichCellsActive[j];
            _concs[iCurr][cell] = conc;
        }

        for (int j = 0; j < dirichCellsActive.size(); j++) {
            auto cell = dirichCellsActive[j];
            freeVector[cell] = conc * _local->_alphas[cell];
        }
    }

}

void Equation::calcConcsImplicit() {

    BiCGSTAB biCGSTAB;

    biCGSTAB.compute(matrix);

    _concs[iCurr] = biCGSTAB.solveWithGuess(freeVector, _concs[iPrev]);

}

void Equation::calcConcsExplicit() {}

void Equation::cfdProcedureOneStep(const double &timeStep) {

    std::swap(iCurr, iPrev);
    _convective->calcBetas(_concs[iPrev]);
    _local->calcAlphas(_concs[iPrev], timeStep);

    processNonBoundFaces(_sgrid->_typesFaces.at("active_nonbound"));
    fillMatrix();
    processDirichCells(_boundGroupsDirich, _concsBoundDirich);

    calcConcsImplicit();

}

void Equation::cfdProcedure() {

    _concs.emplace_back(_sgrid->_cellsArrays.at("concs_array1"));
    _concs.emplace_back(_sgrid->_cellsArrays.at("concs_array2"));

    _local->calcTimeSteps();

    for (auto &timeStep : _local->_timeSteps) {
        cfdProcedureOneStep(timeStep);
        Eigen::Map<Eigen::VectorXd> concCurr(new double[dim], dim);
        concCurr = _concs[iCurr];
        _concsTime.push_back(concCurr);
    }
}

double Equation::calcFacesFlowRate(Eigen::Ref<Eigen::VectorXui64> faces) {

    auto &isMatrix = _sgrid->_cellsConditions.at("is_matrices");

    auto poroFrac = std::get<double>(_props->_params["poro_frac"]);
    auto poroMatrix = std::get<double>(_props->_params["poro_matrix"]);

    auto dFreeFrac = std::get<std::vector<double>>(_props->_params["d_free_frac"]);
    auto dFreeMatrix = std::get<std::vector<double>>(_props->_params["d_free_matrix"]);

    double totalFlowRate = 0;
    for (uint64_t i = 0; i < faces.size(); i++) {
        auto &face = faces[i];

        auto &neighborsCells = _sgrid->_neighborsCells[face];
        auto &cell0 = neighborsCells[0];
        auto &cell1 = neighborsCells[1];

        auto &conc_prev0 = _concs[iPrev](cell0);
        auto &conc_prev1 = _concs[iPrev](cell1);

        auto &normalsNeighborsCells = _sgrid->_normalsNeighborsCells[face];
        auto &norm0 = normalsNeighborsCells[0];
        auto &norm1 = normalsNeighborsCells[1];

        auto diffusivity0 = calcDFree(conc_prev0, isMatrix[cell0], dFreeFrac, dFreeMatrix);
        auto diffusivity1 = calcDFree(conc_prev1, isMatrix[cell1], dFreeFrac, dFreeMatrix);

        auto &axis = _sgrid->_facesAxes[face];
        auto diffusivity = _convective->weighing("meanAverage",
                                                 diffusivity0, diffusivity1);

        auto &conc_curr0 = _concs[iCurr](cell0);
        auto &conc_curr1 = _concs[iCurr](cell1);

        auto poro0 = calcPoro(conc_prev0, poroFrac, poroMatrix, isMatrix[cell0]);
        auto poro1 = calcPoro(conc_prev1, poroFrac, poroMatrix, isMatrix[cell1]);

        auto poro = _convective->weighing("meanAverage", poro0, poro1);

        auto &dS = _sgrid->_facesSs[axis];
        auto &dL = _sgrid->_spacing[axis];

        totalFlowRate -= diffusivity * poro *
                         (norm0 * conc_curr0 + norm1 * conc_curr1) * dS / dL;
    }

    return totalFlowRate;
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

void Equation::setConcsTime(
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