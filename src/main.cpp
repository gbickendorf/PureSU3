#include "LatticeSettings.h"
#include "Matrix.h"
#include "Random.h"
#include "SUNLattice.h"
#include "Vectors.h"
#include <complex>
#include <float.h>
#include <inttypes.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <sys/resource.h>
#include <utility>
#include <vector>
#define Complex complex<double>
using namespace std;
using namespace SUN;

//Routine to calculate the temperature dependence of correllators
void TempPlot() {
  char name[100];
  int T = 10;
  int L = 10;
  vector<int> autocorr{1,   1,   1,   1,   1,   2,   2,   3,   3,
                       4,   7,   10,  15,  25,  47,  92,  184, 259,
                       333, 401, 460, 512, 567, 614, 644, 674, 701,
                       727, 738, 771, 778, 797, 810, 821};
  vector<string> names;
  vector<LatticeSettings> settings;
  vector<SUNLattice *> lats;

  //Define types of plaquettes investigated
  vector<pair<int, int>> prototype;
  prototype.push_back(pair<int, int>(1, 1));
  prototype.push_back(pair<int, int>(1, 2));
  prototype.push_back(pair<int, int>(1, 3));
  prototype.push_back(pair<int, int>(1, 4));
  prototype.push_back(pair<int, int>(1, 5));
  prototype.push_back(pair<int, int>(2, 2));
  prototype.push_back(pair<int, int>(2, 3));

//Prepare lattices of different temperatures
  for (int i = 0; i < 34; i++) {
    LatticeSettings set((i + 1) / 3.0, L, T, StartCondition::Hot);
    settings.push_back(set);
  }
  for (int i = 0; i < 34; i++) {
    sprintf(name, "DAT/Wilson%i.dat", i);
    lats.push_back(new SUNLattice(settings[i]));
    names.push_back(name);
  }
  //Run MC simulation
#pragma omp parallel for schedule(static, 1)
  for (int i = 0; i < 34; i++) {
    lats[i]->WilsonLoopRun(names[i].c_str(), autocorr[i], 1000, prototype);
    printf("%s----Done\n\n", names[i].c_str());
  }
}

// Routine to study finite size effects of the lattice calculation
void FiniteSize() {

  char name[100];
  vector<int> autocorr{5, 4, 5, 3, 4, 3, 254, 261, 269, 253, 261, 265};
  vector<string> names;
  vector<LatticeSettings> settings;
  vector<SUNLattice *> lats;
  vector<pair<int, int>> prototype;
  prototype.push_back(pair<int, int>(1, 1));
  prototype.push_back(pair<int, int>(1, 2));
  prototype.push_back(pair<int, int>(1, 3));
  prototype.push_back(pair<int, int>(2, 2));
  prototype.push_back(pair<int, int>(2, 3));

  // Create latticeproperties of variing sizes and temperature of 3 and 6
  for (int i = 3; i < 10; i++) {
    settings.push_back(LatticeSettings(3.0, i, i, StartCondition::Hot));
    settings.push_back(LatticeSettings(6.0, i, i, StartCondition::Hot));
  }

  // Create actual lattices
  for (int i = 0; i < 12; i++) {
    sprintf(name, "DAT/FiniteSize/FiniteSize%i.dat", i);
    lats.push_back(new SUNLattice(settings[i]));
    names.push_back(name);
  }
  // Run MonteCarlo simulation on all 12 configurations
  #pragma omp parallel for
  for (int i = 0; i < 12; i++) {
    lats[i]->WilsonLoopRun(names[i].c_str(), autocorr[i], 10000, prototype);
    printf("%s----Done\n\n", names[i].c_str());
  }
}

// Routine to extract behaviour of autocorrelation time
void autocorrPlot() {
  LatticeSettings settings(8.0, 10, 10, StartCondition::Hot);
  SUNLattice lat(settings);
  lat.Run("autocorr10", 10);
}

int main() {
  FiniteSize();
  return 0;
}
