#include "numeric.h"

double Newton::solveEq(std::function<double( double& )> EnergyCalc, double startValue)
{
    double dh = 0.0001;
    double t1 = startValue + dh, t2 = startValue;

    auto foo = EnergyCalc;

    double f1 = foo(t1);
    double f2 = foo(t2);
    double df = (f1 - f2) / (dh);

    double nextT = startValue - foo(startValue) / df; // первое приближение
    double prevT  = startValue; // первое приближение
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

