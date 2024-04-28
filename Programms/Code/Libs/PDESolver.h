//
// Created by Иван on 4/28/2024.
//

#ifndef CODE_PDESOLVER_H
#define CODE_PDESOLVER_H

#include "PDEProblem.h"
#include "algebra.h"
#include "FileIO.h"

bool CrossScheme(const PDEProblem &problem, const std::string &filename ="CrossScheme");

#endif //CODE_PDESOLVER_H
