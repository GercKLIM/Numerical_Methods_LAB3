cmake_minimum_required(VERSION 3.27)
project(Code)

set(CMAKE_CXX_STANDARD 14)


add_executable(Code
        Libs/algebra.cpp
        Libs/algebra.h
        Libs/FileIO.cpp
        Libs/FileIO.h
        InputData
        OutputData
        Libs/PDEProblem.cpp
        Libs/PDEProblem.h
        Libs/PDESolver.cpp
        Libs/PDESolver.h
        main.cpp
)
