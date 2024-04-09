#pragma once
#include "macroparam.h"
#include "energycalc.h"
#include <functional>

struct NonLinearEqSolver
{
    virtual double solveEq(std::function<double( double& )> EnergyCalc, double startValue) = 0;
};
struct Newton : public NonLinearEqSolver
{
    double solveEq(std::function<double( double& )> EnergyCalc, double startValue);
};
