#include "LatticeSettings.h"
#include <complex>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define Complex complex<double>
using namespace std;

#ifndef VECTORS_H
#define VECTORS_H

//Define properies a vector on the lattice needs.
class LatticePosition {
public:
  int Pos[4];
  LatticeSettings &settings;
  LatticePosition(int t, int x, int y, int z, LatticeSettings &Settings);

private:
  int wrapt(int t);
  int wraps(int s);
};

LatticePosition operator*(LatticePosition &, int n);
LatticePosition operator+(LatticePosition p1, LatticePosition p2);
LatticePosition operator-(LatticePosition p1, LatticePosition p2);
#endif
