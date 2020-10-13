#include "funcs.h"

double calcAFunc(const double &conc, const double &poro) {
    return poro;
}

double calcBFunc(const double &conc, const double &poro) {
    double diffusivity = 0.002;
    return poro * diffusivity;
}

