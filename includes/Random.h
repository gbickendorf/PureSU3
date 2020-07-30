#include "LatticeSettings.h"
#include "Matrix.h"
#include "Vectors.h"
#include <complex>
#include <cstdint>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef RANDOM_H
#define RANDOM_H

#define Complex complex<double>
using namespace std;
//Define properties the RNG needs for MC lattice calculations
class Random {
private:
  LatticeSettings &settings;
  mt19937_64 mt;
  uniform_real_distribution<double> dist_real01;

public:
  Random(LatticeSettings &Settings);
  double GetDouble();
  Complex GetComplex();
  int GetInt(int limit);
  matrix2 GetUnitaryMatrix2();
  matrix3 GetUnitaryMatrix3();
  LatticePosition GetLatticePosition();
};
;

#endif
