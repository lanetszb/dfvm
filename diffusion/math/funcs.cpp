#include "funcs.h"

double calcDFree(const bool &isMatrix,
                 const std::vector<double> &dFreeFrac,
                 const std::vector<double> &dFreeMatrix) {

    if (isMatrix)
        return dFreeMatrix[1];
    else return dFreeFrac[1];

}

double calcDSurface(const bool &isMatrix,
                    const std::vector<double> &dSurfaceMatrix) {

    if (isMatrix)
        return dSurfaceMatrix[1];
    else return dSurfaceMatrix[0];

}

double calcPoro(const double &conc, const double &poroFrac, const double &poroMatrix,
                const bool &isMatrix) {

    if (isMatrix)
        return poroMatrix;
    else return poroFrac;

}

double calcAFunc(const double &conc, const double &poroFrac, const double &poroMatrix,
                 const bool &isMatrix) {

    return calcPoro(conc, poroFrac, poroMatrix, isMatrix);
}

double calcBFunc(const double &conc, const double &poroFrac, const double &poroMatrix,
                 const bool &isMatrix,
                 const std::vector<double> &dFreeFrac,
                 const std::vector<double> &dFreeMatrix,
                 const std::vector<double> &dSurfaceMatrix) {

    double b;

    if (isMatrix)
        b = calcDFree(isMatrix, dFreeFrac, dFreeMatrix) + calcDSurface(isMatrix, dSurfaceMatrix);
    else
        b = calcDFree(isMatrix, dFreeFrac, dFreeMatrix);

    return b;
}

