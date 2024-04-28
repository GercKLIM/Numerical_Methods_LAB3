//
// Created by Иван on 4/28/2024.
//

#ifndef CODE_PDEPROBLEM_H
#define CODE_PDEPROBLEM_H

#include <cmath>
#include <iostream>
#include <functional>

class PDEProblem {
public:
    double rho; // плотность материала
    double k;   // коэффициент натяжения
    double a;   // скорость распространения малых возмущений
    double x0;  // начало отсчёта по пространству (координата левого конца)
    double L;   // характерный пространственный рамзер (длина струны)
    double X;   // координата правого конца
    double t0;  // время отсчёта
    double T;   // время окончания
    double tau; // шаг по времени
    double h;   // шаг по пространству
    double gam; // число Куранта
    int num_time_steps; // количество шагов по времени
    int num_space_steps; // количество шагов по пространству
    PDEProblem(double rho_init, double k_init, double x0_init, double L_init, double t0_init, double T_init, double h_init, double tau_init);
    PDEProblem(double a_init, double x0_init, double L_init, double t0_init, double T_init, double h_init, double tau_init, bool what_is_L_init);
    // Отклонение точки в нулевой момент времени (f(x))
    std::function<double(double)> initDeflectionFunc;

    // Скорость точки в нулевой момент времени (g(x))
    std::function<double(double)> initVelocityFunc;

    // Закон движения левого конца струны (\varphi(t))
    std::function<double(double)> leftBoundaryFunction;

    // Закон движения правого конца струны (\psi(t))
    std::function<double(double)> rightBoundaryFunction;

    // Воздействие внешних сил (F(x,t)) (если == 0, то уравнение - однородное)
    std::function<double(double, double)> extForcesFunction;
};


#endif //CODE_PDEPROBLEM_H
