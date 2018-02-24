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
	vector<string> names;
	vector<LatticeSettings> settings;
	vector<SUNLattice*> lats;
	names.push_back("DAT/wilson1.dat");
	names.push_back("DAT/wilson2.dat");
	names.push_back("DAT/wilson3.dat");
	names.push_back("DAT/wilson4.dat");
	names.push_back("DAT/wilson5.dat");
	names.push_back("DAT/wilson6.dat");
	names.push_back("DAT/wilson7.dat");
	names.push_back("DAT/wilson8.dat");
	settings.push_back(LatticeSettings(1.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(2.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(3.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(4.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(5.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(6.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(7.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(8.0,6,6,StartCondition::Hot));
	//SUNLattice lat = SUNLattice(settings[0]);
	for (int i = 0; i < 8; i++)
	{
		lats.push_back(new SUNLattice(settings[i]));
	}
	#pragma omp parallel for
	for (int i = 0; i < 8; i++)
	{
		lats[i]->Run(names[i].c_str());
	}
	
	
	//lats.push_back(SUNLattice( LatticeSettings(3.0,6,6,StartCondition::Hot)));
	//LatticeSettings set = LatticeSettings(3.0,6,6,StartCondition::Hot);
	//SUNLattice lat(set);
	//lat.Run("DAT/wilson1.dat");
	return 0;
}
