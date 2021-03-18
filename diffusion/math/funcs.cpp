#include "funcs.h"

double calcPoro(const double &conc, const double &poro, const bool &isMatrix) {

    if (isMatrix)
        return poro;
    else return poro;
}

double calcAFunc(const double &conc, const double &poro, const bool &isMatrix) {

    double a;

    if (isMatrix)
        a = poro;
    else
        a = poro;

    return a;
}

double calcBFunc(const double &conc, const double &dFreeFrac, const double &dFreeMatrix,
                 const double &dSurfaceMatrix, const double &poro, const bool &isMatrix) {

    double b;

    if (isMatrix)
        b = dFreeMatrix + dSurfaceMatrix;
    else
        b = dFreeFrac;

    return b;
}

