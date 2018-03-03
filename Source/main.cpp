#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include <sys/resource.h>
#include <vector>
#include <utility>
#include <iostream>
#include <complex>
#include "Matrix.h"
#include "Random.h"
#include "Vectors.h"
#include "LatticeSettings.h"
#include "SUNLattice.h"
#define Complex complex<double>
using namespace std;
using namespace SUN;


int main()
{
	char name[100];
	int T = 16;
	int L = 16;
	vector<int> autocorr{1,1,2,2,3,5,10,20,45,95,154,203};
	vector<string> names;
	vector<LatticeSettings> settings;
	vector<SUNLattice*> lats;
	vector<pair<int,int>> prototype;
	prototype.push_back(pair<int,int>(1,1));
	prototype.push_back(pair<int,int>(1,2));
	prototype.push_back(pair<int,int>(2,2));
	prototype.push_back(pair<int,int>(2,3));
	prototype.push_back(pair<int,int>(3,3));
	prototype.push_back(pair<int,int>(3,4));
	prototype.push_back(pair<int,int>(4,4));
	settings.push_back(LatticeSettings(1.0,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(1.5,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(2.0,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(2.5,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(3.0,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(3.5,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(4.0,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(4.5,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(5.0,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(5.4,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(5.7,L,T,StartCondition::Hot));
	settings.push_back(LatticeSettings(6.5,L,T,StartCondition::Hot));
	for (int i = 0; i < 12; i++)
	{
		
		sprintf(name,"DAT/Wilson%i.dat",i);
		lats.push_back(new SUNLattice(settings[i]));
		names.push_back(name);
	}
	#pragma omp parallel for
	for (int i = 0; i < 12; i++)
	{
		lats[i]->WilsonLoopRun(names[i].c_str(),autocorr[i], 10000, prototype);
	}
	return 0;
}
