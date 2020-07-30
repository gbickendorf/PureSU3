#include "Matrix.h"

matrix2::matrix2(MatrixInitType init) {
  (*this).M[0][0] = 0.0;
  (*this).M[1][0] = 0.0;
  (*this).M[0][1] = 0.0;
  (*this).M[1][1] = 0.0;
  if (init == MatrixInitType::Unit) {
    (*this).M[1][1] = 1.0;
    (*this).M[0][0] = 1.0;
  }
}

matrix2::matrix2() { matrix2(MatrixInitType::Zero); }

//Prints dim 2 matrix nicely
void matrix2::print() {
  printf("-------------------------------\n");
  printf("| %1.3f%+1.3fi | %1.3f%+1.3fi |\n", (*this).M[0][0].real(),
         (*this).M[0][0].imag(), (*this).M[0][1].real(),
         (*this).M[0][1].imag());
  printf("| %1.3f%+1.3fi | %1.3f%+1.3fi |\n", (*this).M[1][0].real(),
         (*this).M[1][0].imag(), (*this).M[1][1].real(),
         (*this).M[1][1].imag());
  printf("-------------------------------\n");
}

//dagger operation on matrix (complex conjugate)
matrix2 matrix2::dagger() {
  matrix2 dag;
  dag.M[0][0] = conj((*this).M[0][0]);
  dag.M[0][1] = conj((*this).M[1][0]);
  dag.M[1][0] = conj((*this).M[0][1]);
  dag.M[1][1] = conj((*this).M[1][1]);
  return dag;
}

//Calculate trace of dim 2 matrix
Complex matrix2::trace() { return (*this).M[0][0] + (*this).M[1][1]; }

//Calculate determinant of dim 2 matrix
Complex matrix2::det() {
  return (*this).M[0][0] * (*this).M[1][1] - (*this).M[1][0] * (*this).M[0][1];
}

//Defined matrix addition of dim2 matrix
matrix2 operator+(matrix2 m1, matrix2 m2) {
  matrix2 sum;
  sum.M[0][0] = m1.M[0][0] + m2.M[0][0];
  sum.M[0][1] = m1.M[0][1] + m2.M[0][1];
  sum.M[1][0] = m1.M[1][0] + m2.M[1][0];
  sum.M[1][1] = m1.M[1][1] + m2.M[1][1];
  return sum;
}

//Defined matrix subtraction of dim2 matrix
matrix2 operator-(matrix2 m1, matrix2 m2) {
  matrix2 sum;
  sum.M[0][0] = m1.M[0][0] - m2.M[0][0];
  sum.M[0][1] = m1.M[0][1] - m2.M[0][1];
  sum.M[1][0] = m1.M[1][0] - m2.M[1][0];
  sum.M[1][1] = m1.M[1][1] - m2.M[1][1];
  return sum;
}

//Defined matrix multiplication of dim2 matrix
matrix2 operator*(matrix2 m1, matrix2 m2) {
  matrix2 prod;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        prod.M[i][j] += m1.M[i][k] * m2.M[k][j];
      }
    }
  }
  return prod;
}

matrix3::matrix3(MatrixInitType init) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*this).M[i][j] = 0.0;
    }
  }

  if (init == MatrixInitType::Unit) {
    (*this).M[2][2] = 1.0;
    (*this).M[1][1] = 1.0;
    (*this).M[0][0] = 1.0;
  }
}
 //default initialisation type 0
matrix3::matrix3() { matrix3(MatrixInitType::Zero); }

//Prints dim 3 matrix nicely
void matrix3::print() {
  printf("-------------------------------\n");
  printf("| %1.3f%+1.3fi | %1.3f%+1.3fi | %1.3f%+1.3fi |\n",
         (*this).M[0][0].real(), (*this).M[0][0].imag(), (*this).M[0][1].real(),
         (*this).M[0][1].imag(), (*this).M[0][2].real(),
         (*this).M[0][2].imag());
  printf("| %1.3f%+1.3fi | %1.3f%+1.3fi | %1.3f%+1.3fi |\n",
         (*this).M[1][0].real(), (*this).M[1][0].imag(), (*this).M[1][1].real(),
         (*this).M[1][1].imag(), (*this).M[1][2].real(),
         (*this).M[1][2].imag());
  printf("| %1.3f%+1.3fi | %1.3f%+1.3fi | %1.3f%+1.3fi |\n",
         (*this).M[2][0].real(), (*this).M[2][0].imag(), (*this).M[2][1].real(),
         (*this).M[2][1].imag(), (*this).M[2][2].real(),
         (*this).M[2][2].imag());
  printf("-------------------------------\n");
}

//dagger operation on dim 3 matrix
matrix3 matrix3::dagger() {
  matrix3 dag;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dag.M[i][j] = conj((*this).M[j][i]);
    }
  }

  return dag;
}

//Calculate trace of dim 3 matrix
Complex matrix3::trace() {
  return (*this).M[0][0] + (*this).M[1][1] + (*this).M[2][2];
}

//Calculate determinant of dim 3 matrix
Complex matrix3::det() {
  return (*this).M[0][0] * (*this).M[1][1] * (*this).M[2][2] +
         (*this).M[0][1] * (*this).M[1][2] * (*this).M[2][0] +
         (*this).M[0][2] * (*this).M[1][0] * (*this).M[2][1] -
         (*this).M[0][2] * (*this).M[1][1] * (*this).M[2][0] -
         (*this).M[0][1] * (*this).M[1][0] * (*this).M[2][2] -
         (*this).M[0][0] * (*this).M[1][2] * (*this).M[2][1];
}

//Defined matrix addition of dim 3 matrix
matrix3 operator+(matrix3 m1, matrix3 m2) {
  matrix3 sum;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sum.M[i][j] = m1.M[i][j] + m2.M[i][j];
    }
  }
  return sum;
}
//Defined matrix subtraction of dim 3 matrix
matrix3 operator-(matrix3 m1, matrix3 m2) {
  matrix3 sum;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sum.M[i][j] = m1.M[i][j] - m2.M[i][j];
    }
  }
  return sum;
}

//Defined matrix multiplication of dim 3 matrices
matrix3 operator*(matrix3 m1, matrix3 m2) {
  matrix3 prod;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        prod.M[i][j] += m1.M[i][k] * m2.M[k][j];
      }
    }
  }
  return prod;
}
//Defined scalar-matrix multiplication of dim 3 matrix
matrix3 operator*(Complex c, matrix3 m2) {
  matrix3 prod;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      prod.M[i][j] = c * m2.M[i][j];
    }
  }
  return prod;
}
