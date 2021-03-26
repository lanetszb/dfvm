#ifndef DFVM_FUNCS_H
#define DFVM_FUNCS_H

#include <iostream>
#include <map>
#include <vector>

#include <Eigen/Dense>

double calcConcSurfacePrime(const double &conc, const bool &isMatrix);

double calcDFree(const double &conc, const bool &isMatrix,
                 const std::vector<double> &dFreeFrac,
                 const std::vector<double> &dFreeMatrix);

double calcDSurface(const double &conc, const bool &isMatrix,
                    const std::vector<double> &dSurfaceMatrix);

double calcPoro(const double &conc, const double &poroFrac, const double &poroMatrix,
                const bool &isMatrix);

double calcAFunc(const double &conc, const double &poroFrac, const double &poroMatrix,
                 const bool &isMatrix);

double calcBFunc(const double &conc, const double &poroFrac, const double &poroMatrix,
                 const bool &isMatrix,
                 const std::vector<double> &dFreeFrac,
                 const std::vector<double> &dFreeMatrix,
                 const std::vector<double> &dSurfaceMatrix);

#endif //DFVM_FUNCS_H
