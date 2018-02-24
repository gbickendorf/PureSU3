#include "Random.h"

vector<std::mt19937_64> mt_;
vector<std::uniform_real_distribution<double>> dist_real01_;

Random::Random(LatticeSettings & Settings):settings(Settings)
{
	mt = mt19937_64(12);
	dist_real01 = uniform_real_distribution<double>(0,1);
}

double Random::GetDouble()
{
	return dist_real01(mt);
}

int Random::GetInt(int limit)
{
	return uniform_int_distribution<int>{0,limit-1}(mt);
}

Complex Random::GetComplex()
{
	return Complex(GetDouble()-0.5,GetDouble()-0.5);
}

matrix2 Random::GetUnitaryMatrix2()
{
	matrix2 mat;
	Complex c1 = GetComplex();
	Complex c2 = GetComplex();
	double magn = sqrt(norm(c1)+norm(c2));
	c1/=magn;
	c2/=magn;
	mat.M[0][0]=c1;
	mat.M[1][1]=conj(c1);
	mat.M[0][1]=c2;
	mat.M[1][0]=-conj(c2);
	if(GetDouble() <0.5)
		return mat;
	else
		return mat.dagger();
}

matrix3 Random::GetUnitaryMatrix3()
{
	Complex c1[3] ={GetComplex(),GetComplex(),GetComplex()};
	Complex c2[3] ={GetComplex(),GetComplex(),GetComplex()};
	double n1=sqrt(norm(c1[0])+norm(c1[1])+norm(c1[2]));
	
	c1[0]/=n1;
	c1[1]/=n1;
	c1[2]/=n1;
	
	Complex n2 =c2[0]*conj(c1[0])+c2[1]*conj(c1[1])+c2[2]*conj(c1[2]);
	c2[0]-=n2*c1[0];
	c2[1]-=n2*c1[1];
	c2[2]-=n2*c1[2];
	n2=sqrt(norm(c2[0])+norm(c2[1])+norm(c2[2]));
	
	c2[0]/=n2;
	c2[1]/=n2;
	c2[2]/=n2;
	
	matrix3 mat;
	mat.M[0][0]=c1[0];
	mat.M[0][1]=c1[1];
	mat.M[0][2]=c1[2];
	mat.M[1][0]=c2[0];
	mat.M[1][1]=c2[1];
	mat.M[1][2]=c2[2];
	mat.M[2][0]=conj(c1[1]*c2[2]-c1[2]*c2[1]);
	mat.M[2][1]=conj(c1[2]*c2[0]-c1[0]*c2[2]);
	mat.M[2][2]=conj(c1[0]*c2[1]-c1[1]*c2[0]);
	if(GetDouble()<0.5)
		return mat;
	else
		return mat.dagger();
	
}

LatticePosition Random::GetLatticePosition()
{
	return LatticePosition(GetInt(settings.t_LattSize),GetInt(settings.s_LattSize),GetInt(settings.s_LattSize),GetInt(settings.s_LattSize),settings);
}
