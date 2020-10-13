#ifndef DFVM_FUNCS_H
#define DFVM_FUNCS_H

#include <iostream>
#include <map>
#include <vector>

#include <Eigen/Dense>

double calcAFunc(const double &conc, const double &poro);

double
calcBFunc(const double &conc, const double &diffusion, const double &poro);


#endif //DFVM_FUNCS_H
