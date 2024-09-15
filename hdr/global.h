/**
 * @file global.h
 * @brief Definitions of global constants and some basic structures
 * 
 * This file contains the definitions of various global constants used throughout the solver,
 * as well as utility structures and classes, such as solverParameters structure and Mixture class, which
 * is a and vector class wrapper.
 */

#pragma once

#include <vector>
#include <cmath>
#include <string>

using namespace std;

// Main constants

const double N_A = 6.02214076E23; /**< Avogadro's number (1/mol) */
const double R_U = 8.3144598;     /**< Universal Gas constant (J/mol/K) */
const double kB = 1.38064852e-23; /**< Boltzmann's constant (J/K) */
const double kBE = 8.617e-5;      /**< Boltzmann's constant (elVolt*K) */

static const double hPlank = 6.62607015e-34; /**< Planck's constant (J*s) */
static const double clight = 2.99792458e8;   /**< Speed of light (m/s) */
static const double hc = hPlank * clight; 

/**
 * @brief Parameters for the solver
 * 
 * This structure holds various parameters used in the solver.
 */
struct solverParams
{ 
    int NumCell = 0;          /**< Number of computational cells */
    double Ma = 0;            /**< Mach number */
    double Gamma = 0;         /**< Adiabatic index */
    double CFL = 0;           /**< Courant number */
    double lambda = 0;        /**< Mean free path length */
    int MaxIter = 10000000;   /**< Maximum number of time steps */
};

/**
 * @brief Additional wrapper class on vector 
 * 
 * This class is a wrapper on std::vector<double> with additional methods. Currently, it is not used in the solver.
 */

/*
class Matrix
{
private:
	std::vector<double> data;
public:
	Matrix(vector<double> val)
	{
		data = val;
	}
	Matrix()
	{

	}
	void clear()
	{
		data.clear();
	}
	Matrix(int len, double val = 0)
	{
		data = vector<double>(len, val);
	}
	vector<double>::iterator begin()
	{
		return data.begin();
	}
	vector<double>::iterator end()
	{
		return data.end();
	}
	int size()
	{
		return data.size();
	}
	operator vector<double>() const
	{
		return data;
	}
	double& operator [](int i)
	{
		static double elseValue;
		if (i < data.size())
			return data[i];
		else
			return elseValue;
	}
	double first()
	{
		return data.front(); 
	}
	double last()
	{
		return data.back();
	}
	void push_back(double val)
	{
		data.push_back(val);
	}
	void push_front(double val)
	{
		data.insert(data.begin(), val); 
	}
	void removeLast()
	{
		data.pop_back();
	}
	void removeFirst()
	{
		data.erase(data.begin());
	}
	void resize(int i)
	{
		data.resize(i);
	}
	Matrix operator /(vector<double> div)
	{
		vector<double> output;
		if (data.size() != div.size())
			return data;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] / div[i]);
		return output;
	}
	Matrix operator *(vector<double> div)
	{
		vector<double> output;
		if (data.size() != div.size())
			return data;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] * div[i]);
		return output;
	}
	Matrix operator +(vector<double> div)
	{
		vector<double> output;
		if (data.size() != div.size())
			return data;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] + div[i]);
		return output;
	}
	Matrix operator -(vector<double> div)
	{
		vector<double> output;
		if (data.size() != div.size())
			return data;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] - div[i]);
		return output;
	}
	Matrix operator /(double div)
	{
		vector<double> output;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] / div);
		return output;
	}
	Matrix operator *(double div)
	{
		vector<double> output;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] * div);
		return output;
	}
	Matrix operator +(double div)
	{
		vector<double> output;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] + div);
		return output;
	}
	Matrix operator -(double div)
	{
		vector<double> output;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] - div);
		return output;
	}
	static Matrix POW(vector<double> div, double param)
	{
		vector<double> output;
		for (int i = 0; i < div.size(); i++)
			output.push_back(pow(div[i], param));
		return output;
	}
	static Matrix SQRT(vector<double> div)
	{
		vector<double> output;
		for (int i = 0; i < div.size(); i++)
			output.push_back(sqrt(div[i]));
		return output;
	}
	static Matrix REVERSE(vector<double> div)
	{
		vector<double> output;
		for (int i = 0; i < div.size(); i++)
			output.push_back(1.0 / div[i]);
		return output;
	}
};
*/
