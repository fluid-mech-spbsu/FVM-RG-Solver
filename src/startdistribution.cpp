/**
 * @file startdistribution.cpp
 * @brief Implementation of the StartDistribution and its derived classes functions.
 * 
 * This file contains the implementation of setStartDistribution functions that set initial distributions of macroscopic parameters.
 */

#include"startdistribution.h"


void UniformDistribution::setStartDistribution(vector<macroParam>& points)
{
	points[0].mixture = example.mixture;
	points[points.size() - 1].mixture = example.mixture;
	for (size_t i = 1; i < points.size() - 1; i++)
	{
		points[i].mixture = example.mixture;
		points[i].temp = example.temp;
		points[i].fractionArray = example.fractionArray;
		points[i].density = example.density;

		points[i].pressure = points[i].density * R_U * example.temp / example.mixture.molarMass(example.fractionArray);
		points[i].densityArray = example.densityArray;
		points[i].velocity_tau = example.velocity_tau;
		points[i].velocity_normal = example.velocity_normal;
		points[i].velocity = fabs(points[i].velocity_tau);
	}
	
	borderCondition->updatePoints(points);
	return;
}

void GapDistribution::setStartDistribution(vector<macroParam>& points)
{
	Mixture mixture = left.mixture;
	bool startType = 1;
	if (startType) {
		for (size_t i = 0; i < points.size() / 2 + 1; i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = left.pressure;
			points[i].density = left.density;
			points[i].fractionArray = left.fractionArray;
			points[i].densityArray = left.densityArray;
			points[i].velocity_tau = left.velocity_tau;
			points[i].velocity_normal = left.velocity_normal;
			points[i].velocity = left.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass(left.fractionArray) / points[i].density / R_U;
			points[i].gamma = energyCalculator->getGamma(points[i]);
		}
		for (size_t i = points.size() / 2 + 1; i < points.size(); i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = right.pressure;
			points[i].density = right.density;
			points[i].fractionArray = right.fractionArray;
			points[i].densityArray = right.densityArray;
			points[i].velocity_tau = right.velocity_tau;
			points[i].velocity_normal = right.velocity_normal;
			points[i].velocity = right.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass(right.fractionArray) / points[i].density / R_U;
			points[i].gamma = energyCalculator->getGamma(points[i]);
		}
	}
	else {
		for (size_t i = 0; i < 1; i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = left.pressure;
			points[i].density = left.density;
			points[i].fractionArray = left.fractionArray;
			points[i].densityArray = left.densityArray;
			points[i].velocity_tau = left.velocity_tau;
			points[i].velocity_normal = left.velocity_normal;
			points[i].velocity = left.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass(left.fractionArray) / points[i].density / R_U;
			points[i].gamma = energyCalculator->getGamma(points[i]);
		}
		for (size_t i = 1; i < points.size(); i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = right.pressure;
			points[i].density = right.density;
			points[i].fractionArray = right.fractionArray;
			points[i].densityArray = right.densityArray;
			points[i].velocity_tau = right.velocity_tau;
			points[i].velocity_normal = right.velocity_normal;
			points[i].velocity = right.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass(right.fractionArray) / points[i].density / R_U;
			points[i].gamma = energyCalculator->getGamma(points[i]);
		}
	}
}

void FixedDistribution::setStartDistribution(vector<macroParam>& points)
{
	points = exampleVec;
	return;
}