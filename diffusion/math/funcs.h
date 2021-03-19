#ifndef DFVM_FUNCS_H
#define DFVM_FUNCS_H

#include <iostream>
#include <map>
#include <vector>

#include <Eigen/Dense>

double calcDSurface(const bool &isMatrix,
                    const std::vector<double> &dSurfaceMatrix);

double calcDFree(const bool &isMatrix,
                 const std::vector<double> &dFreeFrac,
                 const std::vector<double> &dFreeMatrix);

double calcPoro(const double &conc, const double &poro, const bool &isMatrix);

double calcAFunc(const double &conc, const double &poro, const bool &isMatrix);

double calcBFunc(const double &conc, const double &poro, const bool &isMatrix,
                 const std::vector<double> &dFreeFrac,
                 const std::vector<double> &dFreeMatrix,
                 const std::vector<double> &dSurfaceMatrix);


#endif //DFVM_FUNCS_H
