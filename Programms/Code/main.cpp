#include <iostream>
#include "Libs/PDESolver.h"

int main()
{
    double a = 1.;
    double x0 = 0.;
    double L = 1.;
    bool what_is_L = true; // флаг на указание размера (в случае false - будем передавать координату правого конца)
    double t0 = 0.;
    double T = 1.;
    double h = 0.02;
    double tau = 0.001;

    // Тест 1
    PDEProblem problem1(a, x0, L, t0, T, h, tau, what_is_L);
    problem1.initDeflectionFunc = ([](double x) {return sin(M_PI * x);});
    problem1.initVelocityFunc = ([](double x){return 0.;});
    problem1.leftBoundaryFunction = ([](double t){return 0.;});
    problem1.rightBoundaryFunction = ([](double t){return 0.;});
    problem1.f_xx = ([](double x) {return -M_PI*M_PI*sin(M_PI * x);});
    problem1.f_xx_is_set = true;

    CrossScheme(problem1, "test1.txt");

    // Тест 2
    PDEProblem problem2(a, x0, L, t0, T, h, tau, what_is_L);
    problem2.initDeflectionFunc = ([](double x) {return x*(1.-x);});
    problem2.initVelocityFunc = ([](double x){return 0.;});
    problem2.leftBoundaryFunction = ([](double t){return 0.;});
    problem2.rightBoundaryFunction = ([](double t){return 0.;});
    problem2.f_xx = ([](double x) {return -2.;});
    problem2.f_xx_is_set = true;

    CrossScheme(problem2, "test2.txt");
    return 0;
}
