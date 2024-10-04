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
    int caseType = 1;
    double T_up_wall;
    double T_down_wall;
    double velocity_up;
    double velocity_down;
    if(caseType == 0)
    {
        T_up_wall = 273;
        T_down_wall = 273;
        velocity_up = 300;
        velocity_down = 0;
    }
    else if(caseType == 1)
    {
        T_up_wall = 273;
        T_down_wall = 273;
        velocity_up = 1888.84;
        velocity_down = 0;
    }


    

    //////////////////////////////////////////////////////////////
    ///////////////////// Boundary Conditions for Couette ///////////
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

    UniformDistribution startParamCouetteAr;
    
    macroParam startParamAr(Ar);
    
    startParamAr.density = 0.00012786; // 0.03168;
    startParamAr.fractionArray[0] = 1;
    startParamAr.densityArray[0] =  startParamAr.fractionArray[0] * startParamAr.density;

    startParamAr.temp = 273; //140
    startParamAr.velocity_tau = 5;
    startParamAr.velocity_normal = 0;

    startParamCouetteAr.setBorderCondition(&borderConditionCouette);
    startParamCouetteAr.setDistributionParameter(startParamAr);

    //////////////////////////////////////////////////////////////

    solverParams solParam;
    solParam.NumCell     = 67;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.67;   // Ar
    solParam.CFL      = 0.9;    // Число Куранта 0.9
    solParam.MaxIter     = 10000000; // максимальное кол-во итареций
    solParam.Ma       = 0.1;    // Число маха

    double precision = 1E-6; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(10000);


    // DataWriter writer(outputData);
    DataWriter writer(outputData);
    DataReader reader(outputData + "/prev-data");

    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);


    double viscocity_argon = 2.2724e-05;
    double pressure = startParamAr.density * R_U * startParamAr.temp / argon.molarMass;
    double MFP = viscocity_argon / pressure * sqrt(M_PI * R_U * startParamAr.temp / argon.molarMass); // mean free path length for argon
    std::cout << "mean free path: " << MFP << std::endl;

    GodunovSolver solver(Ar, solParam, SystemOfEquationType::couette2, RiemannSolverType::HLLESolver);
    double h = 1;

    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));

    solver.setBorderConditions(&borderConditionCouette);
    solver.setStartDistribution(&startParamCouetteAr);
    
    solver.solve();
}