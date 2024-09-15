/**
 * @file mixture.cpp
 * @brief Implementation of the Mixture class functions.
 * 
 * This file contains the implementation of various functions for the Mixture class, 
 * including calculations for molar mass, mass, and other properties of the mixture components.
 */

#include "mixture.h"
#include "global.h"

double Mixture::molarMass(std::vector<double> y_c)
{
	double sum = 0;
	for (size_t j = 0; j < y_c.size(); j++)
	{
		double M_c = molarMass(j);
		sum += y_c[j] / M_c;
	}
	return 1 / sum;
}

double Mixture::molarMass(size_t i)
{
	return components[i].molarMass;
}

double Mixture::mass(size_t i)
{
	return components[i].mass;
}

double Mixture::sigma(size_t i)
{
	return components[i].sigma;
}

double Mixture::epsilonDevK(size_t i)
{
	return components[i].epsilonDevK;
}

double Mixture::Zinf(size_t i)
{
	return components[i].Zinf;
}