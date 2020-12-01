#include "funcs.h"

double calcPoro(const double &conc, const double &poro) {

    return poro;
}

double calcAFunc(const double &conc, const double &poro) {

    return calcPoro(conc, poro);
}

double calcBFunc(const double &conc, const double &diffusivity, const double &poro) {
    return poro * diffusivity;
}

