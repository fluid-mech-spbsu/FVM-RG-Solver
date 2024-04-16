#include "numeric.h"

double Newton::solveEq(std::function<double( double& )> EnergyCalc, double startValue)
{
    double dh = 0.000001;
    double t1 = startValue + dh, t2 = startValue;

    auto foo = EnergyCalc;

    double f1 = foo(t1);
    double f2 = foo(t2);
    double df = (f1 - f2) / (dh);

    double nextT = startValue - foo(startValue) / df; // первое приближение
    double prevT  = startValue;
    double eps = dh/100;
    while (fabs(nextT - prevT) > eps)
    {
        prevT = nextT;
        t1 = prevT + dh;
        t2 = prevT;
        f1 = foo(t1);
        f2 = foo(t2);
        df = (f1 - f2) / (dh);

        nextT = prevT - foo(prevT) / df; // последующие приближения
    }
    return nextT;
}


double IntegralCalculatorTrap::integrate(std::vector<double> x, std::vector<double> y)
{
    double sum = 0;
    for(size_t i = 1; i < y.size(); i++)
    {
        sum += (y[i-1] + y[i]) / 2. * (x[i] - x[i-1]);
    }
    return sum;
}
double IntegralCalculatorTrap::integrate(std::vector<double> y, double dh)
{
    double sum = 0;
    for(size_t i = 1; i < y.size(); i++)
    {
        sum += (y[i-1] + y[i]) / 2. * dh;
    }
    return sum;
}

double DerrivativeCalculator::derrivativeCenter(std::vector<double> y, double dh, int idx)
{
    if(idx == 0 || idx == y.size() - 1)
        return 0;
    return (y[idx + 1] - y[idx - 1]) / (2 * dh);
}

double DerrivativeCalculator::derrivativeLeft(std::vector<double> y, double dh, int idx)
{
    if(idx == 0)
        return 0;
    return (y[idx] - y[idx - 1]) / (dh);
}

double DerrivativeCalculator::derrivativeRight(std::vector<double> y, double dh, int idx)
{
    if(idx == y.size() - 1)
        return 0;
    return (y[idx + 1] - y[idx]) / (dh);
}
