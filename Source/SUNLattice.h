#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <complex>
#include <iostream>
#include "LatticeSettings.h"
#include "Matrix.h"
#include "Vectors.h"
#include "Random.h"
#include "LatticeSettings.h"
#include "TwoLinkOperator.h"
#define Complex complex<double>

#define forlatt for (int z = 0; z < settings.s_LattSize; z++)\
				for (int y = 0; y < settings.s_LattSize; y++)\
				for (int x = 0; x < settings.s_LattSize; x++)\
				for (int t = 0; t < settings.t_LattSize; t++)
#define GROUPDIM 3.0
using namespace std;

#ifndef SUNLATT_H
#define SUNLATT_H
namespace SUN{
class SUNLattice
{
public:
	SUNLattice(LatticeSettings& Settings);
	~SUNLattice();
	void Run(const char * filename, int autocorrTime);
	int CoordinatesToIndex(int t, int x, int y, int z, int mu);
	int LatticePositionToIndex(LatticePosition pos, int mu);
	LatticePosition IndexToLatticePosition(int n, int * mu);
	void WilsonLoopRun(const char * filename,int autocorrTime, int N, const vector<pair<int,int>> Loops);
	LatticeSettings& settings;
private:
	int acc, dec;
	Random* random;
	matrix3* Lattice;
	matrix3 GetLink(LatticePosition pos, int mu);
	matrix3 Plaquette(LatticePosition n, int mu, int nu, int a, int b);
	matrix3	Staple(LatticePosition n, int mu);
	double LocalSChange(LatticePosition n, int mu, matrix3 &);
	void MonteCarloStep();
	void MonteCarloUpdate();
	double WilsonLoop(int a, int b);
	void GenQuarkPot();
	Complex PlaqPlaqCorrelation(int tdiff);
	Complex PolyakovLoop(int x, int y, int z);
	void PolyakovLoopCorrelation();
	void PolyakovRun();
	void PlaqPlaqRun();
	int PolyakovPosToIndex(int x, int y, int z);
	void PolyakovIndexToPos(int index, int &, int &, int &);
	void AvgActionPerPlaqueRun();
	void SaveLattice();
	void LoadLattice();	
	int AutoCorrelationWilson(int N, int a, int b);
	void Thermalize(int Ntherm, int verbose);
	
//MultiAlg
	int ThreeCoordinatesToIndex(int x, int y, int z, int timeSlice, int direction, int SliceCount);
	void AvgLevel(int level);
	int ThreeCoordinatesToIndex(int x, int y, int z, int direction);
	void Avg1aSlices(int slice, TwoLinkOperator * result);
	void Avg2aSlices(int slice, TwoLinkOperator * result);
	void Avg4aSlices(int slice, TwoLinkOperator * result);
	double TraceAvgComlete();
	void MultiRun();
	void MultiRun2();
	void TimeSiliceUpdate(int x0, int y0);
	Complex TwoLinkOp(LatticePosition x, int spaceDir, int seperation, int a,int b, int c ,int d);
};
}
#endif
