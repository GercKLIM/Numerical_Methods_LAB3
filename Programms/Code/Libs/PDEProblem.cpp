//
// Created by Иван on 4/28/2024.
//


#include "PDEProblem.h"

PDEProblem::PDEProblem(double rho_init, double k_init, double x0_init, double L_init, double t0_init, double T_init, double h_init, double tau_init) {
    rho = rho_init;
    k = k_init;
    if(rho != 0)
        a = k/rho;
    else{
        throw std::invalid_argument( "received zero rho!" );
    }
    x0 = x0_init;
    L = L_init;
    t0 = t0_init;
    T = T_init;
    h = h_init;
    tau = tau_init;
    gam = fabs(a * tau / h);
    if(gam > 1)
        std::cout << "LOG[WARNING]: Courant number is greater than 1!!!" << std::endl;
    X = x0 + L;
    num_time_steps = static_cast<int>((T-t0) / tau);
    num_space_steps = static_cast<int>((X - x0)/h);
}

PDEProblem::PDEProblem(double a_init, double x0_init, double L_init, double t0_init, double T_init, double h_init, double tau_init, bool what_is_L_init = true) {
    x0 = x0_init;
    if(what_is_L_init){
        L = L_init;
        X = x0 + L;
    }
    else {
        X = L_init;
        L = X-x0;
    }
    t0 = t0_init;
    T = T_init;
    h = h_init;
    tau = tau_init;
    gam = fabs(a * tau / h);
    if(gam > 1)
        std::cout << "LOG[WARNING]: Courant number is greater than 1!!!" << std::endl;

    num_time_steps = static_cast<int>((T-t0) / tau);
    num_space_steps = static_cast<int>((X - x0)/h);
    // Мнимые значения для rho и k (для инициализации)
    rho = -1.;
    k = -1.;
}
