#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <iostream>
#include <cstdint>
#include <random>
#include <omp.h>
#include "Matrix.h"
#include "Vectors.h"
#include "LatticeSettings.h"

#ifndef RANDOM_H
#define RANDOM_H

#define Complex complex<double>
using namespace std;

class Random
{
private: 
	LatticeSettings& settings;
	mt19937_64 mt;
	uniform_real_distribution<double> dist_real01;
public: 
	Random(LatticeSettings& Settings);
	double GetDouble();
	Complex GetComplex();
	int GetInt(int limit);
	matrix2 GetUnitaryMatrix2();
	matrix3 GetUnitaryMatrix3();
	LatticePosition GetLatticePosition();
};
;

#endif
