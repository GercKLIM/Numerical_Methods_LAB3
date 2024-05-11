#include <iostream>
#include "Libs/PDESolver.h"


void konevs_test() {

    double eps = 1e-10;
    // Сетка по пространству и времени (h, tau)
    std::vector<vector<double>> Web = {
          /*  {0.1, 0.01}, // gam = 0.1
            {0.05, 0.025}, // gam = 0.5
            {0.025, 0.01875},// gam = 0.75
            {0.01, 0.01} // gam = 1.0*/
             {0.1, 0.01}, // gam = 0.1
            {0.1, 0.05}, // gam = 0.5
            {0.1, 0.075},// gam = 0.75
            {0.1, 0.1}   // gam = 1.0
    };

    for (int i = 0; i < 4; i++) {
        /* Задача 1 */
        double a = 1.;
        double x0 = -2.;
        double L = 2.;
        bool what_is_L = false;
        double t0 = 0.;
        double T = 1.;
        double h = Web[i][0];
        double tau = Web[i][1];


        PDEProblem problem1(a, x0, L, t0, T, h, tau, what_is_L);

        // f(x)
        problem1.initDeflectionFunc = ([](double x) {                   // f(x)
            if ((x < 1. / 3) and (x > -1. / 3)) {
                return 1;
            } else {
                return 0;
            }
        });

        problem1.initVelocityFunc = ([](double x) { return 0.; });       // g(x)
        problem1.leftBoundaryFunction = ([](double t) { return 0.; });   // phi(x)
        problem1.rightBoundaryFunction = ([](double t) { return 0.; });  // psi(x)
        //problem1.f_xx = ([&](double x) { return 0; });                      // f''_xx анал

        // f''_xx числ
        problem1.f_xx = ([&](double x) { return (problem1.initDeflectionFunc(x + eps) + 2 * problem1.initDeflectionFunc(x) + problem1.initDeflectionFunc(x - eps)) / (eps * eps);});
        problem1.f_xx_is_set = true;
        //problem1.reflective_boundaries = true;

        CrossScheme(problem1, "konev_test1_" + to_string(problem1.gam) + ".txt");


        /* Задача 2*/
        a = 1.;
        x0 = -1.;
        L = 1.;
        what_is_L = false;
        t0 = 0.;
        T = 1.;

        h = Web[i][0];
        tau = Web[i][1];

        PDEProblem problem2(a, x0, L, t0, T, h, tau, what_is_L);

        problem2.initDeflectionFunc = ([](double x) { return 0; });       // f(x)

        problem2.initVelocityFunc = ([](double x) {                    // g(x)
            if ((x > -0.5) and (x < 0.5)) {
                return (1. - fabs(x));
            } else {
                return 0.;
            };
        });
        problem2.leftBoundaryFunction = ([](double t) { return 0.; });   // phi(x)
        problem2.rightBoundaryFunction = ([](double t) { return 0.; });  // psi(x)
        //problem2.f_xx = ([&](double x) { return 0.; });                 // f''_xx анал

        // f''_xx числ
        problem2.f_xx = ([&](double x) {return (problem2.initDeflectionFunc(x + eps) + 2 * problem2.initDeflectionFunc(x) + problem2.initDeflectionFunc(x - eps)) / (eps * eps);});
        problem2.f_xx_is_set = true;
        //problem2.reflective_boundaries = true;

        CrossScheme(problem2, "konev_test2_" + to_string(problem2.gam) + ".txt");

        /* Задача 3*/
        a = 1.;
        x0 = 0.;
        L = 4. * M_PI;
        what_is_L = true;
        t0 = 0.;
        T = 1.;

        h = Web[i][0];
        tau = Web[i][1];

        PDEProblem problem3(a, x0, L, t0, T, h, tau, what_is_L);

        problem3.initDeflectionFunc = ([](double x) { return 0.; });   // f(x)

        problem3.initVelocityFunc = ([](double x) { return 0.; });
        problem3.leftBoundaryFunction = ([](double t) { return sin(t); });   // phi(x)
        problem3.rightBoundaryFunction = ([](double t) { return 0.; });  // psi(x)
        //problem3.f_xx = ([&](double x) { return 0.; });                 // f''_xx анал

        // f''_xx числ
        problem3.f_xx = ([&](double x) {return (problem3.initDeflectionFunc(x + eps) + 2 * problem3.initDeflectionFunc(x) + problem3.initDeflectionFunc(x - eps)) / (eps * eps); });
        problem3.f_xx_is_set = true;

        CrossScheme(problem3, "konev_test3_" + to_string(problem3.gam) + ".txt");
    }
}

/* Функция для создания данных для таблиц */
void make_data_for_tables() {

    /* Тестируем на тесте №1 из методички */
    double a = 1.;
    double x0 = 0.;
    double L = 1.;
    bool what_is_L = true; // флаг на указание размера (в случае false - будем передавать координату правого конца)
    double t0 = 0.;
    double T = 1.;
    double h = 0.2;
    double tau = 0.1;

    // Тест 1


    for (int i = 0; i < 6; i++) {

        PDEProblem problem1(a, x0, L, t0, T, h, tau, what_is_L);
        problem1.initDeflectionFunc = ([](double x) {return sin(M_PI * x);});
        problem1.initVelocityFunc = ([](double x){return 0.;});
        problem1.leftBoundaryFunction = ([](double t){return 0.;});
        problem1.rightBoundaryFunction = ([](double t){return 0.;});
        problem1.f_xx = ([](double x) {return -M_PI*M_PI*sin(M_PI * x);});
        problem1.f_xx_is_set = true;

        tau = tau / 4;
        h = h / 4;
        CrossScheme(problem1, "data_for_tables/test1/test1_" + to_string(i) + ".txt");
    }
}

/* Тесты из методички */
void tests() {
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
}
int main() {
    //tests();
    make_data_for_tables();
    //konevs_test();
    return 0;

}
