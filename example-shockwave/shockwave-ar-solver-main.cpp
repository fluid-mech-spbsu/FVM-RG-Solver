#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>
#include <energycalc.h>
#include <chrono>


std::string GetCurrentWorkingDir(void) {
	std::filesystem::path currentWorkingDir = std::filesystem::current_path();
	std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

	std::string res = parentDir.string() + "/FVM-RG-Solver/example-shockwave";
	return res;
}

namespace fs = std::filesystem;
int main()
{

	std::string outputData = GetCurrentWorkingDir();
	std::cout << "Current directory is: " << outputData << std::endl;

	///////////////////////////////////////////////////////////////
	///////////////////// Particle data //////////////////////////
	///////////////////////////  Ar  ////////////////////////////
	
	MixtureComponent argon;
    argon.name = "Ar";
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21/kB; 
    argon.sigma = 3.33E-10;
    argon.numberAtoms = 1;
	argon.numberOfModes = 0;
	argon.omega_e = 0;

	std::vector<MixtureComponent> tmp = {argon};
    Mixture Ar(tmp);

	// double Pr = 0.664; // Prandtl number for argon
	// Mixture Ar(tmp, Pr);

	//////////////////////////////////////////////////////////////
	//////////////////// Approximation, 1T ///////////////////////
	//////////////////////////////////////////////////////////////

	OneTempApprox oneTempApprox;

	//////////////////////////////////////////////////////////////
	///////////////// Boundary Conditions for Shock Wave /////////
	////////////////////////// CH4 ///////////////////////////////

	double velocity_left = 1225.88; // m/s, Mach 3.8
	double density_left =  0.0016015; // kg/m^3, when p = 100 Pa
	double T_left = 300; // K
	double pressure_left = R_U * T_left * density_left / argon.molarMass;
    double viscocity_argon = 2.2724e-05; // Pa*s, for the above set of conditions

    // Calculate mean free path length for Argon
    double MFP = viscocity_argon / pressure_left * sqrt(M_PI * R_U * T_left / argon.molarMass); // Mean free path length for Argon
    std::cout << "Mean free path: " << MFP << std::endl;

	// General relations, from python solver
	// double velocity_right =  296.11053; 
	// double density_right = 0.00663031; 
	// double T_right = 1395.252; 

	// Rankine-Hugoniot boundary conditions, from python solver 
	double velocity_right = 370.248;
	double density_right = 0.005302675;
	double T_right = 1612.935; 

	double pressure_right = R_U * T_right * density_right / argon.molarMass;

	BorderConditionShockwave borderConditionShockwave;
	borderConditionShockwave.setBorderParameters(
		velocity_left, density_left, T_left,
		velocity_right, density_right, T_right
	);

	borderConditionShockwave.setEnergyCalculator(&oneTempApprox);

	//////////////////////////////////////////////////////////////
	////////////////// Start params for Shockwave ////////////////
	////////////////////////////  Ar  ////////////////////////////

	ShockwaveGapDistribution startParamShockwaveAr; // The area is divided into two parts, in each of which the gas parameters are constant
	macroParam leftStartParam(Ar);
	macroParam rightStartParam(Ar);

	leftStartParam.density = density_left;
	leftStartParam.pressure = pressure_left;
	leftStartParam.velocity_tau = velocity_left;
	leftStartParam.velocity_normal = 0;
	leftStartParam.velocity = velocity_left;
	leftStartParam.fractionArray[0] = 1;
	leftStartParam.densityArray[0] = leftStartParam.fractionArray[0] * leftStartParam.density;

	rightStartParam.density = density_right;
	rightStartParam.pressure = pressure_right;
	rightStartParam.velocity_normal = 0;
	rightStartParam.velocity_tau = velocity_right;
	rightStartParam.velocity = velocity_right;
	rightStartParam.fractionArray[0] = 1;
	rightStartParam.densityArray[0] = rightStartParam.fractionArray[0] * rightStartParam.density;

	startParamShockwaveAr.setDistributionParameter(leftStartParam, rightStartParam);
	startParamShockwaveAr.setEnergyCalculator(&oneTempApprox);

	//////////////////////////////////////////////////////////////
	///////////////// Solver params for Shockwave ////////////////
	////////////////////////////////////////////////////////////

	solverParams solParam;
    solParam.NumCell = 30 + 2;  // Number of computational cells including two fictitious cells
    solParam.Gamma = 1.67;      // Adiabatic index for Argon
    solParam.CFL = 0.9;         // Courant number
    solParam.MaxIter = 20000;    // Maximum number of iterations
    solParam.Ma = 3.8;          // Mach number (currently not affecting the solver, just a formality)

    // Set precision for the observer
    double precision = 1E-6;    // Precision for the observer
    Observer watcher(precision);
    watcher.setPeriodicity(1000); // Set periodicity for the observer

    // Initialize data writer and reader
    DataWriter writer(outputData);
    DataReader reader(outputData + "/prev_data");

    // Read previous data
    reader.read();
    std::vector<macroParam> startParameters;
    reader.getPoints(startParameters);

    // Initialize the solver
    GodunovSolver solver(Ar, solParam, SystemOfEquationType::shockwave1, RiemannSolverType::HLLESolver);
    solver.setOutputDirectory(outputData);

    std::cout << "mean free path: " << MFP << std::endl;
	double h = 20 * MFP; // m
	std::cout << "considering h = MFP * " << h / MFP << std::endl;
	writer.setDelta_h(h / (solParam.NumCell - 2));
	solver.setWriter(&writer);
	solver.setObserver(&watcher);
	solver.setDelta_h(h / (solParam.NumCell - 2));

	solver.setEnergyCalculator(&oneTempApprox);
	solver.setBorderConditions(&borderConditionShockwave);
	solver.setStartDistribution(&startParamShockwaveAr);

	std::cout << "Start solving\n";
	auto start = std::chrono::high_resolution_clock::now();
	solver.solve();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	double seconds = elapsed.count();
	printf("Final time of solution: %f seconds\n", seconds);

}
