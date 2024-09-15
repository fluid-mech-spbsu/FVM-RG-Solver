/**
 * @file coeffsolver.h
 * @brief Definitions for the CoeffSolver and its derived classes
 * 
 * This file contains the definitions of the CoeffSolver class and its derived classes,
 * which provide various methods for calculating coefficients such as viscosity, thermal conductivity, and diffusion.
 */

#pragma once
#include "global.h"
#include "macroparam.h"
#include "energycalc.h"

/**
 * @brief Abstract base class for coefficient solvers
 * 
 * This class defines the interface for various coefficient calculations.
 */
struct CoeffSolver
{
public:
	virtual double shearViscosityHS(macroParam currentPoint) = 0;
	virtual double shearViscosityLJ(Mixture mix, double currentT) = 0;
	virtual double thermalCond(macroParam currentPoint) = 0;
	virtual double thermalCondConstPr(macroParam currentPoint) = 0;
	virtual double thermalCondMultiAtom(macroParam currentPoint) = 0;
	virtual double bulkViscosity(macroParam currentPoint) = 0;
	virtual double bulkViscosityMultiAtom(macroParam currentPoint) = 0;
	virtual double bulkViscosityConstViscRel(macroParam currentPoint) = 0;

	virtual double getOmega11(Mixture mix, double T) = 0;
	virtual double getOmega22(Mixture mix, double T) = 0;
	
	virtual double effDiffusion(macroParam currentPoint, size_t component) { return 0; };
	virtual double binaryDiffusion(macroParam currentPoint, size_t comp1, size_t comp2) { return 0; };
};

/**
 * @brief Coefficient solver for a single component gas, 1T model
 * 
 * This class implements the coefficient calculations for a single component gas and one-temperature model.
 */
struct CoeffSolver1Comp1Temp : public CoeffSolver
{
	double shearViscosityHS(macroParam currentPoint);
	double shearViscosityLJ(Mixture mix, double currentT);
	double thermalCond(macroParam currentPoint);
	double thermalCondConstPr(macroParam currentPoint);
	double thermalCondMultiAtom(macroParam currentPoint);
	double bulkViscosity(macroParam currentPoint);
	double bulkViscosityMultiAtom(macroParam currentPoint);
	double bulkViscosityConstViscRel(macroParam currentPoint);

	double getOmega11(Mixture mix, double T);
	double getOmega22(Mixture mix, double T);

private:
	OneTempApproxMultiModes OneTempApprox;
};

/**
 * @brief Coefficient solver for a two-component gas, 1T model
 * 
 * This class implements the coefficient calculations for a two-component gas and one-temperature model.
 */
struct CoeffSolver2Comp1Temp : public CoeffSolver1Comp1Temp
{
	double shearViscosityHS(macroParam currentPoint);
	double thermalCond(macroParam currentPoint);
	double thermalCondConstPr(macroParam currentPoint) { return 0; };
	double bulkViscosity(macroParam currentPoint);
	double bulkViscosityMultiAtom(macroParam currentPoint) { return 0; };
	double bulkViscosityConstViscRel(macroParam currentPoint) { return 0; };
	double effDiffusion(macroParam currentPoint, size_t component);
	double binaryDiffusion(macroParam currentPoint, size_t comp1, size_t comp2);
private:
	double Xc(macroParam currentPoint, size_t component);
	double shearViscosity(macroParam currentPoint, size_t component);
	double lambda(macroParam currentPoint, size_t component);
	double phi(macroParam currentPoint, size_t component1, size_t component2);
	double getOmega11(macroParam currentPoint, size_t comp1, size_t comp2);
};
