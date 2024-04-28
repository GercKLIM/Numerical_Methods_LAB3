#include <iostream>
#include "Libs/PDESolver.h"

int main()
{
    double a = 1.;
    double x0 = -2.;
    double L = 4.;
    bool what_is_L = true; // флаг на указание размера (в случае false - будем передавать координату правого конца)
    double t0 = 0.;
    double T = 2.;
    double h = 0.5;
    double tau = 0.05;

    // Тест 1
    PDEProblem problem1(a, x0, L, t0, T, h, tau, what_is_L);
    problem1.initDeflectionFunc = ([](double x) {return sin(M_PI * x);});
    problem1.initVelocityFunc = ([](double x){return 0.;});
    problem1.leftBoundaryFunction = ([](double t){return 0.;});
    problem1.rightBoundaryFunction = ([](double t){return 0.;});

    CrossScheme(problem1, "test1");

    return 0;
}
