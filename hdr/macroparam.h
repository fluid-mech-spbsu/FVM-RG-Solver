/**
 * @file macroparam.h
 * @brief Definition of the macroParam structure
 * 
 * This file contains the definition of the macroParam structure, which is used to store the 
 * macroscopic parameters of the gas mixture.
 */
#pragma once
#include "mixture.h"

/**
 * @brief Structure to store the macroscopic parameters of the gas mixture
 * 
 * This structure is used to store the macroscopic parameters of the gas mixture, 
 * such as density, pressure, velocity, temperature, etc.
 */
struct macroParam
{
	macroParam() { mixture = Mixture(); }
	macroParam(Mixture mix) :mixture(mix) { fractionArray.resize(mix.NumberOfComponents); 
											densityArray.resize(mix.NumberOfComponents); }
	Mixture mixture; /**< Mixture of gases */
	vector<double> densityArray; /**< Array of densities of the components of the mixture */
	vector<double> fractionArray; /**< Array of mass fractions of the components of the mixture */
	double density = 0; /**< Density of the gas */
	double pressure = 0; /**< Pressure of the gas */
	double velocity_tau = 0; /**< Tangential velocity component */
	double velocity_normal = 0; /**< Normal velocity component */
	double velocity = 0; /**< Total velocity */
	double temp = 0; /**< Gas temperature */
	double temp_vibr = 0; /**< Vibrational temperature, currently is not applied */
	double soundSpeed = 0; /**< Speed of sound in the gas mixture */
	double gamma; /**< Specific heat ratio */
	string gas = "Ar"; /**< Gas mixture name */
};
