#include "funcs.h"

double calcPoro(const double &conc, const double &poro, const bool &isMatrix) {

    if (isMatrix)
        return poro;
    else return poro * 0.2;
}

double calcAFunc(const double &conc, const double &poro, const bool &isMatrix) {

    if (isMatrix)
        return calcPoro(conc, poro, isMatrix);
    else calcPoro(conc, poro, isMatrix) * 0.01;
}

double calcBFunc(const double &conc, const double &diffusivity,
                 const double &poro, const bool &isMatrix) {

    if (isMatrix)
        return calcPoro(conc, poro, isMatrix) * diffusivity;
    else return calcPoro(conc, poro, isMatrix) * diffusivity * 0.01;
}

