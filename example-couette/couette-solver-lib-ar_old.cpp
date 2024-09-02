//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>


std::string GetCurrentWorkingDir( void ) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

   std::string res = parentDir.string() + "/FVM-RG-Solver/example-couette";

    return res;
}

namespace fs = std::filesystem;
int main()
{

    std::string outputData = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << outputData << std::endl;

    //////////////////////////////////////////////////////////////
    ///////////////////// Border Condition for Couette ///////////
    //////////////////////////////////////////////////////////////

    double T_up_wall;
    double T_down_wall;
    double velocity_up;
    double velocity_down;

    T_up_wall = 1000;
    T_down_wall = 1000;
    velocity_up = 300;
    velocity_down = 0;

    //////////////////////////////////////////////////////////////
    ///////////////////// Border Condition for Couette ///////////
    /////////////////////////// Slip /////////////////////////////

    BorderConditionCouette borderConditionCouette;
    borderConditionCouette.setWallParameters(velocity_up, velocity_down, T_up_wall, T_down_wall);

    //////////////////////////////////////////////////////////////

    // Ar
    MixtureComponent argon;
    argon.name = "Ar";
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21/kB;
    argon.numberAtoms = 1;
    argon.sigma = 3.33E-10;


    std::vector<MixtureComponent> tmp = {argon};
    Mixture Ar(tmp);

    //////////////////////////////////////////////////////////////
    ///////////////////// Start param for Couette ////////////////
    ////////////////////////////  Ar  ///////////////////////////

    UniformDistributionBorder startParamCouetteAr;

    startParamCouetteAr.setMixture(Ar); 
    macroParam startParamAr(Ar);
    bool newSolving = true;
    if(newSolving)
    {
        startParamAr.density =  0.000478; //0.00000478;
        startParamAr.fractionArray[0] = 1;
        startParamAr.densityArray[0] =  startParamAr.fractionArray[0] * startParamAr.density;

        startParamAr.temp = 800; //140
        startParamAr.velocity_tau = 30;
        startParamAr.velocity_normal = 0;

        startParamCouetteAr.setBorderCondition(&borderConditionCouette);
        startParamCouetteAr.setDistributionParameter(startParamAr);
    }
    else
    {
        DataWriter writer(outputData);
        DataReader reader(outputData + "/prev_data");

        reader.read();
        vector<macroParam> startParameters;
        reader.getPoints(startParameters);

        startParamCouetteAr.setBorderCondition(&borderConditionCouette);
        startParamCouetteAr.setDistributionParameter(startParameters);
    }

    //////////////////////////////////////////////////////////////

    solverParams solParam;
    solParam.NumCell     = 62;    // Число расчетных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.67;   // Ar
    solParam.CFL      = 0.9;    // Число Куранта 0.9
    solParam.MaxIter     = 10000000; // максимальное кол-во итераций
    solParam.Ma       = 0.1;    // Число маха

    double precision = 1E-5; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(10000);


    DataWriter writer(outputData);
    DataReader reader(outputData + "\prev_data");

    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);

    GodunovSolver solver(Ar, solParam, SystemOfEquationType::couette2Alt, RiemannSolverType::HLLESolver);
    // GodunovSolver solver(Ar, solParam, SystemOfEquationType::couette1, RiemannSolverType::HLLESolver);

    double h = 1;
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));

    solver.setBorderConditions(&borderConditionCouette);
    solver.setStartDistribution(&startParamCouetteAr);
    
    solver.solve();
}