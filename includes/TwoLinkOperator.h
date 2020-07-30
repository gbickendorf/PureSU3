#include "Matrix.h"
#include <complex>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define Complex complex<double>
using namespace std;

#ifndef TWOLINK_H
#define TWOLINK_H

//Define TwoLnkOperators
class TwoLinkOperator {
public:
  Complex trace();
  TwoLinkOperator(matrix3 U1, matrix3 U2);
  TwoLinkOperator();
  void DivideByScalar(double s);
  matrix3 u1;
  matrix3 u2;
};

TwoLinkOperator operator*(TwoLinkOperator m1, TwoLinkOperator m2);
TwoLinkOperator operator+(TwoLinkOperator m1, TwoLinkOperator m2);

#endif
