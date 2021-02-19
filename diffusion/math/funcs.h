#ifndef DFVM_FUNCS_H
#define DFVM_FUNCS_H

#include <iostream>
#include <map>
#include <vector>

#include <Eigen/Dense>

double calcPoro(const double &conc, const double &poro, const bool &isMatrix);

double calcAFunc(const double &conc, const double &poro, const bool &isMatrix);

double calcBFunc(const double &conc, const double &diffusivity,
                 const double &poro, const bool &isMatrix);


#endif //DFVM_FUNCS_H
