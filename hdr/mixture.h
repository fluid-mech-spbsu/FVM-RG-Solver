/**
 * @file mixture.h
 * @brief Definitions for the Mixture and MixtureComponent structures.
 * 
 * This file contains the definitions of the Mixture and MixtureComponent structures,
 * which represent the properties of a mixture and its components, respectively.
 */

#ifndef _MIXTURE
#define _MIXTURE
#include <vector>
#include <string>
#include "global.h"

/**
 * @brief Represents a component (particle) of a mixture
 * 
 * This structure holds various properties of a mixture component. 
 * ! Need to be reorganized.
 */
struct MixtureComponent
{
	double molarMass; /**< Molar mass, set in kg/mol */
    double mass; /**< Mass of a particle, set in kg */
    double sigma; /**< Diameter of a particle, set in m */
    double omega_e; /**< Spectroscopic constants, set in m^-1 */
	double D_diss; /**< Dissociation energy of the molecule, set in m^-1 */
    int numberVibrLvl; /**< Number of vibrational levels */
    int numberAtoms; /**< Number of atoms */
    int numberOfModes; /**< Number of degenerate modes */

	double epsilonDevK; /**< Parameter in the LJ potential, set in K */
    double Zinf; /**< Parameter for Park's model */

    std::vector<double> omega_eByMode; /**< Spectroscopic constants for each degenerate mode */
    std::vector<int> numberVibrLvlByMode; /**< Number of vibrational levels for each degenerate mode */
    std::vector<int> dByMode; /**< Degree of degeneracy for each mode */
    std::vector<std::vector<int>> possibleVibrInds; /**< Set of possible vibrational levels of molecule */

    std::string name; /**< Name of the component */
};

/**
 * @brief Represents a mixture of components
 * 
 * This structure holds various properties of a mixture.
 */
struct Mixture
{
	Mixture() { NumberOfComponents = 0; };
	Mixture(std::vector<MixtureComponent> components_) : components(components_) 
	{ NumberOfComponents = components.size(); 
	  Prandtl = 0.75; // arbitrary value
	};
	Mixture(std::vector<MixtureComponent> components_, double Prandtl_) : components(components_), 
	  Prandtl(Prandtl_) { NumberOfComponents = components.size(); };

	std::vector<MixtureComponent> components; /**< List of mixture components */
    int NumberOfComponents; /**< Number of components in the mixture */
    double Prandtl; /**< Prandtl number */

	double molarMass(std::vector<double> y_c);
	double molarMass(size_t i);
	double mass(size_t i);
	double sigma(size_t i);
	double epsilonDevK(size_t i);
	double Zinf(size_t i);
};
#endif
