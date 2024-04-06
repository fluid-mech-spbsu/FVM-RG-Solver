#include "bordercondition.h"
#include <iostream>
static int counter = 0;
void BorderConditionCouette::updatePoints(vector<macroParam> &points)
{
    bool presEq = 0;
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    if(presEq)
    {
        //0
        points[0].mixture = mixture;
        points[0].pressure = points[1].pressure;
        points[0].fractionArray = points[1].fractionArray;
        points[0].velocity_tau = -points[1].velocity_tau + 2.* down_velocity;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp +  2. * down_temp;
        points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (UniversalGasConstant * points[0].temp);

        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];


        //solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].pressure = points[N-2].pressure;
        points[N-1].fractionArray = points[N-2].fractionArray;
        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* up_velocity;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
        points[N-1].temp = -points[N-2].temp +  2.* up_temp;
        points[N-1].density = points[N-1].pressure * mixture.molarMass(points[N-1].fractionArray) / (UniversalGasConstant * points[N-1].temp);

        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[N-1].fractionArray[i];
    }
    else
    {
        //0
        points[0].mixture = mixture;
        points[0].density = points[1].density;
        points[0].fractionArray = points[1].fractionArray;
        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];

        points[0].velocity_tau = -points[1].velocity_tau + 2.* down_velocity;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp +  2. * down_temp;
        points[0].pressure = points[0].density * UniversalGasConstant * points[0].temp / (mixture.molarMass(points[0].fractionArray));


        //solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].density = points[N-2].density;
        points[N-1].fractionArray = points[N-2].fractionArray;
        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[0].fractionArray[i];

        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* up_velocity;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(fabs(points[N-1].velocity_tau),2) + pow(fabs(points[N-1].velocity_normal),2));
        points[N-1].temp = -points[N-2].temp +  2. * up_temp;
        points[N-1].pressure = points[N-1].density * UniversalGasConstant * points[N-1].temp / (mixture.molarMass(points[N-1].fractionArray));
    }
    return;
}

void BorderConditionPersonal::updatePoints(vector<macroParam> &points)
{

    bool presEq = 1;
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    //0
    points[0].mixture = mixture;

    points[0].pressure = points[1].pressure;

    points[0].density = points[1].density;
    points[0].densityArray = points[1].densityArray;
    points[0].fractionArray = points[1].fractionArray;

    points[0].velocity_tau = points[1].velocity_tau;
    //    points[0].velocity_normal = points[1].velocity_normal;
    points[0].velocity_normal = points[1].velocity_normal;


    points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));

    points[0].temp = points[0].pressure / ( points[0].density * UniversalGasConstant) * mixture.molarMass(points[0].fractionArray);

    // дополнительные рассчитываемые величины

    //solParam.NumCell-1
    points[N-1].mixture = mixture;
    if(presEq)
        points[N-1].pressure = points[N-2].pressure;
    else
        points[N-1].density = points[N-2].density;
    points[N-1].densityArray = points[N-2].densityArray;
    points[N-1].fractionArray = points[N-2].fractionArray;
    points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* up_velocity;
    points[N-1].velocity_normal = -points[N-2].velocity_normal;
    points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
    points[N-1].temp = -points[N-2].temp +  2.* up_temp;
    // дополнительные рассчитываемые величины
    auto y_c = points[N-1].fractionArray;
    if(!presEq)
        points[N-1].pressure = points[N-1].density / mixture.molarMass(y_c) * UniversalGasConstant * points[N-1].temp;
    else
        points[N-1].density = points[N-1].pressure * mixture.molarMass(y_c) / (UniversalGasConstant * points[N-1].temp);
}

void BorderConditionShockwave::updatePoints(vector<macroParam>& points)
{
    bool BCtype = 0;
    size_t N = points.size();
    Mixture mixture = points[1].mixture;

    if (BCtype) {
        points[0].velocity_normal = 0;
        points[0].velocity_tau = left_velocity;
        points[0].velocity = left_velocity;

        points[0].densityArray = points[1].densityArray;
        points[0].fractionArray = points[1].fractionArray;

        points[0].pressure = points[1].pressure;
        points[0].density = points[1].density;
        points[0].temp = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (points[0].density * UniversalGasConstant); // из уравнения состояния ид газа

        points[N-1].velocity_normal = 0;
        points[N-1].velocity_tau = right_velocity;
        points[N-1].velocity = right_velocity;

        points[N-1].densityArray = points[N-2].densityArray;
        points[N-1].fractionArray = points[N-2].fractionArray;

        points[N-1].pressure = points[N-2].pressure;
        points[N-1].density = points[N-2].density;
        points[N-1].temp = points[N-1].pressure * mixture.molarMass(points[N-1].fractionArray) / (points[N-1].density * UniversalGasConstant); // из уравнения состояния ид газа
    }

    // ! another option (check if correct)
    else {
        points[0].velocity_normal = 0;
        points[0].velocity_tau = left_velocity;
        points[0].velocity = points[0].velocity_tau;

        points[0].density = left_density;
        points[0].densityArray[0] = left_density;
        points[0].fractionArray = points[1].fractionArray;
        points[0].temp = left_temp;
        points[0].pressure = points[0].density * UniversalGasConstant * points[0].temp / mixture.molarMass(points[0].fractionArray);

        points[N-1].velocity_normal = 0;
        points[N-1].velocity_tau = right_velocity;
        points[N-1].velocity = points[N-1].velocity_tau;

        points[N-1].density = right_density;
        points[N-1].densityArray[0] = right_density;
        points[N-1].fractionArray = points[N-2].fractionArray;
        points[N-1].temp = right_temp;
        points[N-1].pressure = points[N-1].density * UniversalGasConstant * points[N-1].temp / mixture.molarMass(points[N-1].fractionArray) ;
    }

}

void BorderConditionSoda::updatePoints(vector<macroParam> &points)
{
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    points[0].velocity_normal = 0;
    points[0].velocity_tau = points[1].velocity_tau;
    points[0].velocity = points[0].velocity_tau;

    points[0].densityArray = points[1].densityArray;
    points[0].fractionArray = points[1].fractionArray;

    points[0].pressure = points[1].pressure;
    points[0].density = points[1].density;

    points[N-1].velocity_normal = 0;
    points[0].velocity_tau = points[N-2].velocity_tau;
    points[N-1].velocity = points[N-1].velocity_tau;

    points[N-1].densityArray = points[N-2].densityArray;
    points[N-1].fractionArray = points[N-2].fractionArray;

    points[N-1].pressure = points[N-2].pressure;
    points[N-1].density = points[N-2].density;
}

void BorderConditionCouetteSlip::updatePointsStart(vector<macroParam> &points)
{
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    if(presEq)
    {
        //0
        points[0].mixture = mixture;
        points[0].pressure = points[1].pressure;
        points[0].fractionArray = points[1].fractionArray;
        points[0].velocity_tau = -points[1].velocity_tau + 2.* down_velocity;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp +  2. * down_temp;
        points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (UniversalGasConstant * down_temp);
        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];

        // remember for next iteration
        downLast.temp = interp1(points[0].temp, points[1].temp);
        downLast.fractionArray = points[0].fractionArray;
        downLast.density = points[0].density;
        downLast.densityArray = points[0].densityArray;
        downLast.velocity = down_velocity;
        downLast.velocity_tau = down_velocity;
        downLast.velocity_normal = 0;
        downLast.pressure = points[0].pressure;
        downLast.mixture = mixture;


        //solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].pressure = points[N-2].pressure;
        points[N-1].fractionArray = points[N-2].fractionArray;
        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* up_velocity;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
        points[N-1].temp = -points[N-2].temp +  2.* up_temp;
        points[N-1].density = points[N-1].pressure * mixture.molarMass(points[N-1].fractionArray) / (UniversalGasConstant * up_temp);
        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[N-1].fractionArray[i];

        // remember for next iteration
        upLast.temp = interp1(points[N-1].temp, points[N-2].temp);
        upLast.fractionArray = points[N-1].fractionArray;
        upLast.density = points[N-1].density;
        upLast.densityArray = points[N-1].densityArray;
        upLast.velocity = up_velocity;
        upLast.velocity_tau = up_velocity;
        upLast.velocity_normal = 0;
        upLast.pressure = points[N-1].pressure;
        upLast.mixture = mixture;
    }
    else
    {
        //0
        points[0].mixture = mixture;
        points[0].fractionArray = points[1].fractionArray;
        points[0].density = points[1].density;
        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];

        points[0].velocity_tau = -points[1].velocity_tau + 2.* down_velocity;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp + 2 * down_temp;
        points[0].pressure = points[0].density * UniversalGasConstant * down_temp / mixture.molarMass(points[0].fractionArray) ;

        // remember for next iteration
        downLast.temp = interp1(points[0].temp, points[1].temp);
        downLast.fractionArray = points[0].fractionArray;
        downLast.density = points[0].density;
        downLast.densityArray = points[0].densityArray;
        downLast.velocity = down_velocity;
        downLast.velocity_tau = down_velocity;
        downLast.velocity_normal = 0;
        downLast.pressure = points[0].pressure;
        downLast.mixture = mixture;

        //solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].fractionArray = points[N-2].fractionArray;
        points[N-1].density = points[N-2].density;
        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[N-1].fractionArray[i];

        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* up_velocity;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
        points[N-1].temp = -points[N-2].temp +  2.* up_temp;
        points[N-1].pressure = points[N-1].density * UniversalGasConstant * up_temp / mixture.molarMass(points[N-1].fractionArray) ;

        // remember for next iteration
        upLast.temp = interp1(points[N-1].temp, points[N-2].temp);
        upLast.fractionArray = points[N-1].fractionArray;
        upLast.density = points[N-1].density;
        upLast.densityArray = points[N-1].densityArray;
        upLast.velocity = up_velocity;
        upLast.velocity_tau = up_velocity;
        upLast.velocity_normal = 0;
        upLast.pressure = points[N-1].pressure;
        upLast.mixture = mixture;
    }
}

void BorderConditionCouetteSlip::updatePoints(vector<macroParam> &points)
{
    // NOW CORRECT ONLY FOR SINGLE COMPONENT GAS !!!!!!!!!!!!
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    double velocityHalf;
    double temperatureHalf;
    if(presEq)
    {
        //0
        points[0].mixture = mixture;
        points[0].pressure = points[1].pressure;
        points[0].fractionArray = points[1].fractionArray;

        velocityHalf = calcVelocityHalf(points[1], 0, "down");
        temperatureHalf = calcTempHalf(points[1], 0, velocityHalf, "down");

        points[0].velocity_tau = -points[1].velocity_tau + 2.* velocityHalf;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp + 2 * temperatureHalf;
        points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (UniversalGasConstant * temperatureHalf);
        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];

        // remember for next iteration
        downLast.temp = temperatureHalf;
        downLast.fractionArray = points[0].fractionArray;
        downLast.density = points[0].density;
        downLast.densityArray = points[0].densityArray;
        downLast.velocity = velocityHalf;
        downLast.velocity_tau = velocityHalf;
        downLast.velocity_normal = 0;
        downLast.pressure = points[0].pressure;
        downLast.mixture = mixture;

        // solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].pressure = points[N-2].pressure;
        points[N-1].fractionArray = points[N-2].fractionArray;

        velocityHalf = calcVelocityHalf(points[N-2], 0, "up");
        temperatureHalf = calcTempHalf(points[N-2], 0, velocityHalf, "up");

        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* velocityHalf;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
        points[N-1].temp = -points[N-2].temp +  2.* temperatureHalf;
        points[N-1].density = points[N-1].pressure * mixture.molarMass(points[N-1].fractionArray) / (UniversalGasConstant * temperatureHalf);
        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[N-1].fractionArray[i];

        // remember for next iteration
        upLast.temp = temperatureHalf;
        upLast.fractionArray = points[N-1].fractionArray;
        upLast.density = points[N-1].density;
        upLast.densityArray = points[N-1].densityArray;
        upLast.velocity = velocityHalf;
        upLast.velocity_tau = velocityHalf;
        upLast.velocity_normal = 0;
        upLast.pressure = points[N-1].pressure;
        upLast.mixture = mixture;

    }
    else // does not work
    {
        //0
        points[0].mixture = mixture;
        points[0].fractionArray = points[1].fractionArray;
        points[0].density = points[1].density;
        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];

        velocityHalf = calcVelocityHalf(points[1], 0, "down");
        temperatureHalf = calcTempHalf(points[1], 0, velocityHalf, "down");

        points[0].velocity_tau = -points[1].velocity_tau + 2.* velocityHalf;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp + 2 * temperatureHalf;
        points[0].pressure = points[0].density * UniversalGasConstant * points[0].temp / mixture.molarMass(points[0].fractionArray) ;

        // remember for next iteration
        downLast.temp = temperatureHalf;
        downLast.fractionArray = points[0].fractionArray;
        downLast.density = points[0].density;
        downLast.densityArray = points[0].densityArray;
        downLast.velocity = velocityHalf;
        downLast.velocity_tau = velocityHalf;
        downLast.velocity_normal = 0;
        downLast.pressure = points[0].pressure;
        downLast.mixture = mixture;

        //solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].fractionArray = points[N-2].fractionArray;
        points[N-1].density = points[N-2].density;
        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[N-1].fractionArray[i];

        velocityHalf = calcVelocityHalf(points[N-2], 0, "up");
        temperatureHalf = calcTempHalf(points[N-2], 0, velocityHalf, "up");

        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* velocityHalf;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
        points[N-1].temp = -points[N-2].temp +  2.* temperatureHalf;
        points[N-1].pressure = points[N-1].density * UniversalGasConstant * points[N-1].temp / mixture.molarMass(points[N-1].fractionArray) ;

        // remember for next iteration
        upLast.temp = temperatureHalf;
        upLast.fractionArray = points[N-1].fractionArray;
        upLast.density = points[N-1].density;
        upLast.densityArray = points[N-1].densityArray;
        upLast.velocity = velocityHalf;
        upLast.velocity_tau = velocityHalf;
        upLast.velocity_normal = 0;
        upLast.pressure = points[N-1].pressure;
        upLast.mixture = mixture;
    }
}

macroParam BorderConditionCouetteSlip::getWallParam(string side)
{
    macroParam point;
    if(side == "down")
    {
        point = downLast;
    }
    else if(side == "up")
    {
        point = upLast;
    }
    return point;
}

double BorderConditionCouetteSlip::calcVelocityHalf(macroParam p1, size_t component, string side)
{
    // p0 - ghost cell, p1 - real cell
    int i = component;
    double rhoHalf = p1.densityArray[i];
    double m = p1.mixture.mass(i);
    double M = p1.mixture.molarMass(i);
    macroParam point(p1.mixture);
    double wallVelocity = 0;
    if(side == "down")
    {
        point = downLast;
        wallVelocity = down_velocity;
    }
    else if(side == "up")
    {
        point = upLast;
        wallVelocity = up_velocity;
    }
    double mu = coeffSolver->shareViscositySimple(point); // was forgotten...
    double T_last = point.temp;
    double mult = sqrt(2 * M_PI/ m * kB * T_last) * (2 - sigma) * kB * mu *  M / (sigma * UniversalGasConstant * rhoHalf);
    double numerator = -mult * p1.velocity_tau / delta_h + wallVelocity;
    double denominator = 1 - mult / delta_h;
    return numerator / denominator;
}

double BorderConditionCouetteSlip::calcTempHalf(macroParam p1, size_t component, double velocityHalf, string side)
{
    // p0 - ghost cell, p1 - real cell
    int i = component;
    double rhoHalf = p1.densityArray[i];
    double m = p1.mixture.mass(i);
    double M = p1.mixture.molarMass(i);
    macroParam point(p1.mixture);
    double wallTemperature = 0, wallVelocity = 0;
    if(side == "down")
    {
        point = downLast;
        wallVelocity = down_velocity;
        wallTemperature = down_temp;
        wallVelocity = down_velocity;
    }
    else if(side ==  "up")
    {
        point = upLast;
        wallVelocity = up_velocity;
        wallTemperature = up_temp;
        wallVelocity = up_velocity;
    }
    double T_last = point.temp;
    double mult = (2 - sigma) / (2 * sigma) * sqrt((M_PI * m) / (2 * kB * T_last)) * (M * coeffSolver->lambda(point)) / (UniversalGasConstant * rhoHalf);
    double add = m * pow((velocityHalf - wallVelocity),2) / (4 * kB);
    double numerator = -mult * p1.temp / delta_h + wallTemperature + add ;
    double denominator = 1 - mult / delta_h;
    return numerator / denominator;
}

double BorderConditionCouetteSlip::interp1(double value1, double value2)
{
    return ((value1 + value2) / 2.);
}

