#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <iostream>

#define Complex complex<double>
using namespace std;

#ifndef MATRIX_H
#define MATRIX_H


enum class MatrixInitType { Zero, Unit };
  

class matrix2
{
public:
	Complex det();
	Complex trace();
	void print();
	matrix2();
	matrix2(MatrixInitType init);
	matrix2 dagger();
	Complex M[2][2];
};

matrix2 operator*(matrix2 m1, matrix2 m2);
matrix2 operator+(matrix2 m1, matrix2 m2);
matrix2 operator-(matrix2 m1, matrix2 m2);



class matrix3
{
public:
	Complex det();
	Complex trace();
	void print();
	matrix3();
	matrix3(MatrixInitType init);
	matrix3 dagger();
	Complex M[3][3];
};
matrix3 operator*(Complex c, matrix3 m2);
matrix3 operator*(matrix3 m1, matrix3 m2);
matrix3 operator+(matrix3 m1, matrix3 m2);
matrix3 operator-(matrix3 m1, matrix3 m2);

#endif
