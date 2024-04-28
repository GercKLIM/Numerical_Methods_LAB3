//
// Created by Иван on 4/28/2024.
//
#include "PDESolver.h"

/* Аппроксимация первого слоя (отсчёт слоёв с 0) с порядком O(\tau)
 * args: состояние 0 слоя, объект задачи (содержащий все необходимые параметры)
 * return: состояние 1 слоя
*/
std::vector<double> firstLayerApproximation(const std::vector<double> &state0, const PDEProblem &problem){
    std::vector<double> new_state = state0;
    double x_i = problem.x0;

    new_state[0] = problem.leftBoundaryFunction(x_i);
    for(int i = 1; i < problem.num_space_steps; ++i){
        x_i += problem.h;
        new_state[i] = problem.tau * problem.initVelocityFunc(x_i) + state0[i];
    }

    x_i += problem.h;
    new_state[problem.num_space_steps] = problem.rightBoundaryFunction(x_i);

    return new_state;
}

/* Аппроксимация первого слоя (отсчёт слоёв с 0) с порядком O(\tau^2)
 * args: состояние 0 слоя, объект задачи (содержащий все необходимые параметры), f_xx
 * return: состояние 1 слоя
*/
std::vector<double> firstLayerApproximation(const std::vector<double> &state0, const PDEProblem &problem, std::function<double(double)> f_xx){
    std::vector<double> new_state = state0;
    double x_i = problem.x0;
    double multiplier = problem.a * problem.a * problem.tau * problem.tau / 2;

    new_state[0] = problem.leftBoundaryFunction(x_i);
    for(int i = 1; i < problem.num_space_steps; ++i){
        x_i += problem.h;
        new_state[i] = state0[i] + problem.tau * problem.initVelocityFunc(x_i) + multiplier * f_xx(x_i);
    }
    x_i += problem.h;
    new_state[problem.num_space_steps] = problem.rightBoundaryFunction(x_i);

    return new_state;
}

/* Инициализация начального состояния системы (каждой точке пространства ставим в соотвествие величину начального отклонения)
 * args: problem
 * return: initial state
 * */
std::vector<double> initializeState(const PDEProblem &problem){
    std::vector<double> new_state(problem.num_space_steps+1, 0);
    double x_i = problem.x0;

    for(int i = 0; i <= problem.num_space_steps; ++i){
        new_state[i] = problem.initDeflectionFunc(x_i);
        x_i += problem.h;
    }

    return new_state;
}

bool CrossScheme(const PDEProblem &problem, const string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        // Физические параметры
        double a = problem.a;      // скорость распространения малых возмущений
        double x0 = problem.x0;    // начало отсчёта по пространству
        double L = problem.L;      // характерный пространственный размер (длина струны)
        double X = x0 + L;         // координата праовго края
        double t0 = problem.t0;    // время отсчёта
        double T = problem.T;      // время окончания
        double tau = problem.tau;  // шаг по времени
        double h = problem.h;      // шаг по пространству
        double gam = problem.gam;  // число Куранта

        std::cout << "LOG[INFO]: Courant number gam = " << gam << std::endl;
        if(gam > 1){
            std::cout << "LOG[WARN]: Calculation with gam > 1!" << std::endl;
        }

        // Шаги по времени и пространству
        int num_time_steps = problem.num_time_steps;
        int num_space_steps = problem.num_space_steps;

        // Вычисление нулевого слоя
        std::vector<double> state_j = initializeState(problem);

        // Аппроксимация первого слоя ("jp" = "j+1")
        std::vector<double> state_jp = firstLayerApproximation(state_j, problem);

        // Инициализация второго слоя ("jpp" = "j+1+1")
        std::vector<double> state_jpp = state_jp;

        double t_i = t0;
        double x_i = x0;

        // Запись первых двух слоёв в файл
        writeVectorToFile(fpoints, t_i, state_j);
        t_i += tau;
        writeVectorToFile(fpoints, t_i, state_jp);

        // Эволюция системы во времени
        for(int j = 2; j <= num_time_steps; ++j) {
            t_i += tau;

            // Граничные условия слева
            state_jpp[0] = problem.leftBoundaryFunction(t_i);

            // Граничные условия справа
            state_jpp[num_space_steps] = problem.rightBoundaryFunction(t_i);

            // Обход пространства
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                state_jpp[i] = 2*state_jp[i] - state_j[i] + tau*tau*a*a/(h*h) * (state_jp[i+1] - 2*state_jp[i] + state_jp[i-1]);
            }
            // Запись в файл
            writeVectorToFile(fpoints, t_i, state_jpp);

            // j  jp jpp  -> jp j  jpp -> jp jpp j
            state_j.swap(state_jp);
            state_jpp.swap(state_jp);
        }
        fpoints.close();
        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
};