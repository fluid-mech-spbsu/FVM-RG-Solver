/**
 * @file startdistribution.h
 * @brief Definition of the StartDistribution and its derived structures
 * 
 * This file contains the definition of the StartDistribution structure and its derived structures,
 * which are used to set the initial distribution of macroscopic parameters of the gas mixture.
 */

#pragma once

#include "macroparam.h"
#include "bordercondition.h"
#include "energycalc.h"

/**
 * @brief Base structure for setting the initial distribution of macroscopic parameters
 * 
 * This structure is used as a base class for different types of initial distributions.
 */
struct StartDistribution
{
	/**
     * @brief Pure virtual function to set the initial distribution
     * 
     * This function must be implemented by derived classes to set the initial distribution
     * of macroscopic parameters.
     * 
     * @param points Vector of macroscopic parameters to be initialized
     */
	virtual void setStartDistribution(vector<macroParam>& points) = 0;

	/**
     * @brief Set the energy calculator
     * 
     * This function sets the energy calculator to be used for calculating specific properties.
     * 
     * @param energyCalculator_ Pointer to the energy calculator
     */
	virtual void setEnergyCalculator(EnergyCalc* energyCalculator_) { energyCalculator = energyCalculator_; };

protected:
	EnergyCalc* energyCalculator;
};

/**
 * @brief Structure for setting a uniform initial distribution
 * 
 * This structure sets a uniform initial distribution of macroscopic parameters.
 */
struct UniformDistribution : public StartDistribution
{
	void setStartDistribution(vector<macroParam>& points);
	void setBorderCondition(BorderCondition* borderCondition_) { borderCondition = borderCondition_; };
	void setDistributionParameter(macroParam example_) { example = example_; };
protected:
	BorderCondition* borderCondition;
	macroParam example;
};

/**
 * @brief Structure for setting a gap initial distribution
 * 
 * This structure sets a gap initial distribution of macroscopic parameters, with different
 * parameters on the left and right sides.
 */
struct GapDistribution : public StartDistribution
{
	void setStartDistribution(vector<macroParam>& points);
	void setDistributionParameter(macroParam left_, macroParam right_) { left = left_; right = right_; };
protected:
	macroParam left, right;
};

/**
 * @brief Structure for setting a fixed initial distribution
 * 
 * This structure sets a fixed initial distribution of macroscopic parameters.
 */
struct FixedDistribution : public StartDistribution
{
	void setStartDistribution(vector<macroParam>& points);
	void setDistributionParameter(vector<macroParam>& example_) { exampleVec = example_; };
protected:
	vector<macroParam> exampleVec;
};
