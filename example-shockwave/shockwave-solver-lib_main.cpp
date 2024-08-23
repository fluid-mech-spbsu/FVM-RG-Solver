//#include "hllcsolver.h"
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
	///////////////////// Molecula data //////////////////////////
	///////////////////////////  CH4  ////////////////////////////
	///
	MixtureComponent methane;
	methane.name = "CH4";
	methane.molarMass = 0.016043;
	methane.mass = 2.663732314e-26;
	methane.epsilonDevK = 151.4; // have written in .doc file on drive // epsilon/kB
	methane.Zinf = 89.15;
	methane.sigma = 3.737e-10; // m // another data from the paper rot. relaxation: 3.65e-10
	methane.D_diss = 3668582.3189; // m^-1!, converted from 438.86 kJ/mol
	methane.numberAtoms = 5;
	methane.numberOfModes = 4;
	methane.omega_eByMode = { 302550, 158270, 315680, 136740 }; // m^-1
	methane.numberVibrLvlByMode = { 10, 18, 9, 21 }; // { 4,4,4,4 };
	methane.dByMode = { 1, 2, 3, 3 };

	for (int i1 = 0; i1 < methane.numberVibrLvlByMode[0]; i1++)
	{
		for (int i2 = 0; i2 < methane.numberVibrLvlByMode[1]; i2++)
		{
			for (int i3 = 0; i3 < methane.numberVibrLvlByMode[2]; i3++)
			{
				for (int i4 = 0; i4 < methane.numberVibrLvlByMode[3]; i4++)
				{
					double e = (
						methane.omega_eByMode[0] * (i1 + methane.dByMode[0] / 2.) +
						methane.omega_eByMode[1] * (i2 + methane.dByMode[1] / 2.) +
						methane.omega_eByMode[2] * (i3 + methane.dByMode[2] / 2.) +
						methane.omega_eByMode[3] * (i4 + methane.dByMode[3] / 2.)
						);
					if (e < methane.D_diss) {
						std::vector<int> inds = { i1, i2, i3, i4 };
						methane.possibleVibrInds.push_back(inds);
					}
				}
			}
		}
	}

	// methane.possibleVibrInds = {{0, 0, 0, 0}}; // only ground state

	double viscocity_methane = 1.1123e-05; // for below set sonditions
	// double viscocity_argon = 2.2724e-05; // for below set sonditions


	std::vector<MixtureComponent> tmp3 = { methane };
	double Pr = 0.702; // Prandtl number for methane
	// double Pr = 0.724; // Prandtl number for methane
	Mixture CH4(tmp3, Pr);
	// Mixture CH4(tmp3);

	//////////////////////////////////////////////////////////////
	//////////////////// MODEL ///////////////////////////////////
	//////////////////////////////////////////////////////////////
	OneTempApproxMultiModes oneTempApproxMultiModes;

	///////////////////////////////////////////////////////////////
	///////////////// Border Condition for Shock Wave ////////////
	////////////////////////// CH4 ///////////////////////////////
	///
	// рассматриваем уравнения граничных условий, 0 = left,	n = right.

	// METHANE SET:
	double velocity_left = 1710.228; // Mach 3.8
	double density_left =  0.0006431; // (10 Pa) kg/m^3
	double T_left = 300; // Kelvin
	double pressure_left = UniversalGasConstant * T_left * density_left / methane.molarMass;

	// from python solver (general relations)
	double velocity_right = 263.5241; //Mach 3.8 
	double density_right = 0.004174111; // (10 Pa)
	double T_right =  781.843; 

	// from python solver (Rankine-Hugoniot boundary conditions)
	// double velocity_right = 363.894;
	// double density_right = 0.0030228;
	// double T_right = 1054.0323; 

	// from python solver (general relations, only ground state)
	// double velocity_right = 348.20518;
	// double density_right = 0.00031599; // (10 Pa)
	// double T_right = 976.186; 
	double pressure_right = UniversalGasConstant * T_right * density_right / methane.molarMass;

	BorderConditionShockwave borderConditionShockwave;
	borderConditionShockwave.setBorderParameters(
		velocity_left, density_left, T_left,
		velocity_right, density_right, T_right
	);

	borderConditionShockwave.setEnergyCalculator(&oneTempApproxMultiModes);
	//////////////////////////////////////////////////////////////
	////////////////// Start param for Shockwave /////////////////
	////////////////////////////  CH4  ///////////////////////////

	ShockwaveGapDistribution startParamShockwaveCH4;
	macroParam leftStartParam(CH4);
	macroParam rightStartParam(CH4);

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

	startParamShockwaveCH4.setDistributionParameter(leftStartParam, rightStartParam);
	startParamShockwaveCH4.setEnergyCalculator(&oneTempApproxMultiModes);
	//////////////////////////////////////////////////////////////
	///////////////// Solver param for Shockwave ////////////////
	////////////////////////////////////////////////////////////

	solverParams solParam;
	solParam.NumCell = 60 + 2;  // Число расчетных ячеек с учетом двух фиктивных ячеек
	// solParam.Gamma       = 1.67;        // Ar
	// solParam.Gamma       = 1.32;        // O2_O
	solParam.Gamma = 1.30842;     // CH4, but its also implemented changable in macroparam
	solParam.CFL = 0.9;         // Число Куранта
	solParam.MaxIter = 2000;  // максимальное кол-во итераций
	solParam.Ma = 3.8;			// Число Маха, сейчас не влияет на решатель, просто формальность

	double precision = 1E-6;   // точность
	Observer watcher(precision);
	watcher.setPeriodicity(10000);


	DataWriter writer(outputData);
	DataReader reader(outputData + "/prev_data");

	reader.read();
	vector<macroParam> startParameters;
	reader.getPoints(startParameters);


	// GodunovSolver solver(Ar ,solParam, SystemOfEquationType::shockwave1, RiemannSolverType::HLLESolver);
	GodunovSolver solver(CH4, solParam, SystemOfEquationType::shockwave2, RiemannSolverType::HLLESolver);

	solver.setOutputDirectory(outputData);

	// double MFP = viscocity_argon / pressure_left * sqrt(M_PI * UniversalGasConstant * T_left / argon.molarMass); // mean free path length for argon
	double MFP = viscocity_methane / pressure_left * sqrt(M_PI * UniversalGasConstant * T_left / methane.molarMass); // mean free path length for methane
	std::cout << "mean free path: " << MFP << std::endl;
	double h = 30 * MFP; // m
	std::cout << "considering h = MFP * " << h / MFP << std::endl;
	writer.setDelta_h(h / (solParam.NumCell - 2));
	solver.setWriter(&writer);
	solver.setObserver(&watcher);
	solver.setDelta_h(h / (solParam.NumCell - 2));

	solver.setEnergyCalculator(&oneTempApproxMultiModes);
	solver.setBorderConditions(&borderConditionShockwave);

	// solver.setStartDistribution(&startParamShockwaveAr);
	solver.setStartDistribution(&startParamShockwaveCH4);

	std::cout << "Start solving\n";
	auto start = std::chrono::high_resolution_clock::now();
	// clock_t start = clock();
	solver.solve();
	auto end = std::chrono::high_resolution_clock::now();
	// clock_t end = clock();
	std::chrono::duration<double> elapsed = end - start;
	double seconds = elapsed.count();
	// double seconds = (double)(end - start) / CLOCKS_PER_SEC; // ! The time is calculated incorrectly (e.g. the calculation is performed under 15 min, and the timer shows ~ 9 min) 
	printf("Final time of solution: %f seconds\n", seconds);
}
