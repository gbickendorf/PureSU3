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

void Debug()
{
	LatticeSettings settings(15.0,10,10,StartCondition::Hot);
	SUNLattice lat(settings);
	lat.Run("kjasdj",10);
}

void Diag1()
{
		char name[100];
	int T = 10;
	int L = 10;
	vector<int> autocorr{1,1,1,1,1,2,2,3,3,4,7,10,15,25,47,92,184,259,333,401,460,512,567,614,644,674,701,727,738,771,778,797,810,821};
	vector<string> names;
	vector<LatticeSettings> settings;
	vector<SUNLattice*> lats;
	vector<pair<int,int>> prototype;
	prototype.push_back(pair<int,int>(1,1));
	prototype.push_back(pair<int,int>(1,2));
	prototype.push_back(pair<int,int>(1,3));
	prototype.push_back(pair<int,int>(1,4));
	prototype.push_back(pair<int,int>(1,5));
	prototype.push_back(pair<int,int>(2,2));
	prototype.push_back(pair<int,int>(2,3));

	for (int i = 0; i < 34; i++)
	{
		LatticeSettings set((i+1)/3.0,L,T,StartCondition::Hot);
		settings.push_back(set);
	}
	for (int i = 0; i < 34; i++)
	{
		sprintf(name,"DAT/Wilson%i.dat",i);
		lats.push_back(new SUNLattice(settings[i]));
		names.push_back(name);
	}
	#pragma omp parallel for schedule(static,1)
	for (int i = 0; i < 34; i++)
	{
		//WilsonLoopRun(const char * filename,int autocorrTime, int N, const vector<pair<int,int>> Loops)
		lats[i]->WilsonLoopRun(names[i].c_str(),autocorr[i],1000,prototype);
		printf("%s----Done\n\n",names[i].c_str());
		//lats[i]->WilsonLoopRun(names[i].c_str(),autocorr[i], 10000, prototype);
	}
}

void FiniteSize()
{
	
	char name[100];
	vector<int> autocorr{5,4,5,3,4,3,254,261,269,253,261,265};
	vector<string> names;
	vector<LatticeSettings> settings;
	vector<SUNLattice*> lats;
	vector<pair<int,int>> prototype;
	prototype.push_back(pair<int,int>(1,1));
	prototype.push_back(pair<int,int>(1,2));
	prototype.push_back(pair<int,int>(1,3));
	prototype.push_back(pair<int,int>(2,2));
	prototype.push_back(pair<int,int>(2,3));
	settings.push_back(LatticeSettings(3.0,4,4,StartCondition::Hot));
	settings.push_back(LatticeSettings(3.0,5,5,StartCondition::Hot));
	settings.push_back(LatticeSettings(3.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(3.0,7,7,StartCondition::Hot));
	settings.push_back(LatticeSettings(3.0,8,8,StartCondition::Hot));
	settings.push_back(LatticeSettings(3.0,9,9,StartCondition::Hot));
	settings.push_back(LatticeSettings(6.0,4,4,StartCondition::Hot));
	settings.push_back(LatticeSettings(6.0,5,5,StartCondition::Hot));
	settings.push_back(LatticeSettings(6.0,6,6,StartCondition::Hot));
	settings.push_back(LatticeSettings(6.0,7,7,StartCondition::Hot));
	settings.push_back(LatticeSettings(6.0,8,8,StartCondition::Hot));
	settings.push_back(LatticeSettings(6.0,9,9,StartCondition::Hot));
	for (int i = 0; i < 12; i++)
	{
		sprintf(name,"DAT/FiniteSize/FiniteSize%i.dat",i);
		lats.push_back(new SUNLattice(settings[i]));
		names.push_back(name);
	}
	
	//#pragma omp parallel
	//{
	//	lats[5]->Run(names[5].c_str(),1000);
	//	lats[9]->Run(names[6].c_str(),1000);
	//	lats[10]->Run(names[7].c_str(),1000);
	//	lats[11]->Run(names[8].c_str(),1000);
	//}
	//return ;
	#pragma omp parallel for
	for (int i = 0; i < 12; i++)
	{
		lats[i]->WilsonLoopRun(names[i].c_str(),autocorr[i], 10000, prototype);
		printf("%s----Done\n\n",names[i].c_str());
	}
}


void autocorrPlot()
{
	LatticeSettings settings(8.0,10,10,StartCondition::Hot);
	SUNLattice lat(settings);
	lat.Run("kjasdj",10);
}
int main()
{
	autocorrPlot();
	return 0;
}
