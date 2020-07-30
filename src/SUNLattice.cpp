#include "SUNLattice.h"
#define nIMPROVED

namespace SUN {

matrix3 *Lattice;
int acc(0);
int dec(0);
vector<LatticePosition> Basis;
//Constructor to initialise tools
SUNLattice::SUNLattice(LatticeSettings &Settings) : settings(Settings) {
  Basis.push_back(LatticePosition(1, 0, 0, 0, settings));
  Basis.push_back(LatticePosition(0, 1, 0, 0, settings));
  Basis.push_back(LatticePosition(0, 0, 1, 0, settings));
  Basis.push_back(LatticePosition(0, 0, 0, 1, settings));
  random = new Random(settings);
  Lattice = new matrix3[settings.Volume() * 4];
  if (settings.Condition == StartCondition::Cold) {
    //#pragma omp parallel for
    for (int i = 0; i < settings.Volume() * 4; i++) {
      Lattice[i] = matrix3(MatrixInitType::Unit);
    }
  } else {
    //#pragma omp parallel for
    for (int i = 0; i < settings.Volume() * 4; i++) {
      Lattice[i] = random->GetUnitaryMatrix3();
    }
  }
}
//Translate lattice coordinates to the index
int SUNLattice::CoordinatesToIndex(int t, int x, int y, int z, int mu) {
  return mu + 4 * t + 4 * settings.t_LattSize * x +
         4 * settings.t_LattSize * settings.s_LattSize * y +
         4 * settings.t_LattSize * settings.s_LattSize * settings.s_LattSize *
             z;
}
//Translate lattice position to index n and direction mu
int SUNLattice::LatticePositionToIndex(LatticePosition pos, int mu) {
  return SUNLattice::CoordinatesToIndex(pos.Pos[0], pos.Pos[1], pos.Pos[2],
                                        pos.Pos[3], mu);
}

//Translate lattice index n and direction mu to position
LatticePosition SUNLattice::IndexToLatticePosition(int n, int *mu) {
  int p[4];
  p[3] =
      n / (4 * settings.t_LattSize * settings.s_LattSize * settings.s_LattSize);
  n = n % (4 * settings.t_LattSize * settings.s_LattSize * settings.s_LattSize);
  p[2] = n / (4 * settings.t_LattSize * settings.s_LattSize);
  n = n % (4 * settings.t_LattSize * settings.s_LattSize);
  p[1] = n / (4 * settings.t_LattSize);
  n = n % (4 * settings.t_LattSize);
  p[0] = n / (4);
  n = n % 4;
  *mu = n;
  return LatticePosition(p[0], p[1], p[2], p[3], settings);
}

SUNLattice::~SUNLattice() { delete[] Lattice; }

//Get link at position pos and direction mu
matrix3 SUNLattice::GetLink(LatticePosition pos, int mu) {
  return Lattice[SUNLattice::LatticePositionToIndex(pos, mu)];
}

//Define and calculate geometry of a staple
matrix3 SUNLattice::Staple(LatticePosition n, int mu) {
  matrix3 A = matrix3(MatrixInitType::Zero);
  matrix3 B = matrix3(MatrixInitType::Zero);
  for (int nu = 0; nu < 4; nu++) {
    if (nu == mu)
      continue;
    LatticePosition pos = n + Basis[mu];
    A = A +
        GetLink(n + Basis[mu], nu) * (GetLink(n + Basis[nu], mu).dagger()) *
            (GetLink(n, nu).dagger()) +
        (GetLink(pos - Basis[nu], nu).dagger()) *
            (GetLink(n - Basis[nu], mu).dagger()) * GetLink(n - Basis[nu], nu);
#ifdef IMPROVED
    B = B +
        GetLink(n + Basis[mu], mu) * GetLink(n + Basis[mu] * 2, nu) *
            (GetLink(n + (Basis[mu] + Basis[nu]), mu).dagger()) *
            (GetLink(n + Basis[nu], mu).dagger()) * (GetLink(n, nu).dagger()) +
        GetLink(n + Basis[mu], nu) * (GetLink(n + Basis[nu], mu).dagger()) *
            (GetLink(n + (Basis[nu] - Basis[mu]), mu).dagger()) *
            (GetLink(n - Basis[mu], nu).dagger()) * GetLink(n - Basis[mu], mu) +
        GetLink(n + Basis[mu], nu) * GetLink(n + (Basis[mu] + Basis[nu]), nu) *
            (GetLink(n + Basis[mu] * 2, mu).dagger()) *
            (GetLink(n + Basis[nu], nu).dagger()) * (GetLink(n, nu).dagger()) +
        GetLink(n + Basis[mu], mu) *
            (GetLink(n + (Basis[mu] * 2 - Basis[nu]), nu).dagger()) *
            (GetLink(n + (Basis[mu] - Basis[nu]), mu).dagger()) *
            (GetLink(n - Basis[nu], mu).dagger()) * GetLink(n - Basis[nu], nu) +
        (GetLink(n + (Basis[mu] - Basis[nu]), nu).dagger()) *
            (GetLink(n - Basis[nu], mu).dagger()) *
            (GetLink(n - (Basis[nu] + Basis[mu]), mu).dagger()) *
            GetLink(n - (Basis[mu] + Basis[nu]), nu) *
            GetLink(n - Basis[mu], mu) +
        (GetLink(n + (Basis[mu] - Basis[nu]), nu).dagger()) *
            (GetLink(n + (Basis[mu] - Basis[nu] * 2), nu).dagger()) *
            (GetLink(n + Basis[mu] * 2, mu).dagger()) *
            GetLink(n - Basis[nu] * 2, nu) * GetLink(n - Basis[nu], nu);
#endif
  }
  return A;
}

//Returns the local change of the action given the trial candidate
double SUNLattice::LocalSChange(LatticePosition n, int mu, matrix3 &trial) {
  return real(((GetLink(n, mu) - trial) * Staple(n, mu)).trace()) *
         settings.beta / 3.0;
}

//Update one lattice position by MC
void SUNLattice::MonteCarloStep() {
  LatticePosition pos = random->GetLatticePosition();
  int mu = random->GetInt(4);
  matrix3 trial = random->GetUnitaryMatrix3();
  if (random->GetDouble() < exp(-LocalSChange(pos, mu, trial))) {
    acc++;
    Lattice[LatticePositionToIndex(pos, mu)] = trial;
  } else
    dec++;
}

//Rum on average one MCStep per site
void SUNLattice::MonteCarloUpdate() {
  for (int i = 0; i < settings.Volume() * 4; i++) {
    SUNLattice::MonteCarloStep();
  }
}
//Define and calculate geometry of plaques
matrix3 SUNLattice::Plaquette(LatticePosition n, int mu, int nu, int a, int b) {
  matrix3 xSide(MatrixInitType::Unit);
  matrix3 ySide(MatrixInitType::Unit);
  matrix3 xSideBack(MatrixInitType::Unit);
  matrix3 ySideBack(MatrixInitType::Unit);
  for (int i = 0; i < a; i++) {
    xSide = xSide * GetLink(n + Basis[mu] * i, mu);
    xSideBack =
        xSideBack *
        (GetLink(n + Basis[mu] * (a - i - 1) + Basis[nu] * b, mu).dagger());
  }

  for (int i = 0; i < b; i++) {
    ySide = ySide * GetLink(n + Basis[mu] * a + Basis[nu] * i, nu);
    ySideBack = ySideBack * (GetLink(n + Basis[nu] * (b - i - 1), nu).dagger());
  }

  return xSide * ySide * xSideBack * ySideBack;
}

//Define and calculate geometry of wilsonloops
double SUNLattice::WilsonLoop(int a, int b) {
  double res = 0.0;

  forlatt {
    LatticePosition pos(t, x, y, z, settings);
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = 0; nu < 4; nu++) {
        if (nu == mu)
          continue;
        // count++;
        res += real(Plaquette(pos, mu, nu, a, b).trace());
      }
    }
  }
  return res / (GROUPDIM * settings.Volume() * 12);
}

// Calculate plaque correlation with complex-timedifference tdiff
Complex SUNLattice::PlaqPlaqCorrelation(int tdiff) {
  int size = 1;
  double im, re;
  im = 0.0;
  re = 0.0;
  Complex res = 0.0;
  int count = 0;
#pragma omp parallel for reduction(+ : im, re, count) collapse(4)
  forlatt {
    LatticePosition pos(t, x, y, z, settings);
    for (int mu = 1; mu < 4; mu++) {
      for (int nu = 1; nu < 4; nu++) {
        if (nu == mu)
          continue;
        count++;
        res = Plaquette(pos, mu, nu, size, size).trace() *
              Plaquette(pos + Basis[0] * tdiff, mu, nu, size, size).trace();
        im += imag(res);
        re += real(res);
      }
    }
  }
  return Complex(re, im) / (GROUPDIM * GROUPDIM * count);
}

//Calculate PolyakovLoop at xyz
Complex SUNLattice::PolyakovLoop(int x, int y, int z) {
  matrix3 Loop(MatrixInitType::Unit);
  for (int t = 0; t < settings.t_LattSize; t++) {
    LatticePosition pos(t, x, y, z, settings);
    Loop = Loop * GetLink(pos, 0);
  }
  return Loop.trace();
}

//Translate xyz coordinates to lattice indices
int SUNLattice::PolyakovPosToIndex(int x, int y, int z) {
  return (x % settings.s_LattSize) +
         (y % settings.s_LattSize) * settings.s_LattSize +
         (z % settings.s_LattSize) * settings.s_LattSize * settings.s_LattSize;
}
//Translate lattice indices to xyz coordinates
void SUNLattice::PolyakovIndexToPos(int index, int &x, int &y, int &z) {
  z = index / (settings.s_LattSize * settings.s_LattSize);
  index = index % (settings.s_LattSize * settings.s_LattSize);
  y = index / settings.s_LattSize;
  x = index % settings.s_LattSize;
}

//Calculate correlation of PolyakovLoops in the current configuration
void SUNLattice::PolyakovLoopCorrelation() {
  int avgcound = 0;
  //	int maxSeperation = settings.s_LattSize/2;
  double Result;
  Complex *PolyakovLoops =
      new Complex[settings.s_LattSize * settings.s_LattSize *
                  settings.s_LattSize];
  for (int z = 0; z < settings.s_LattSize; z++) {
    for (int y = 0; y < settings.s_LattSize; y++) {
      for (int x = 0; x < settings.s_LattSize; x++) {
        PolyakovLoops[SUNLattice::PolyakovPosToIndex(x, y, z)] =
            SUNLattice::PolyakovLoop(x, y, z);
      }
    }
  }

  Result = 0.0;
  for (int z = 0; z < settings.s_LattSize; z++) {
    for (int y = 0; y < settings.s_LattSize; y++) {
      for (int x = 0; x < settings.s_LattSize; x++) {
        Result += real(
            PolyakovLoops[SUNLattice::PolyakovPosToIndex(x, y, z)] *
            conj(PolyakovLoops[SUNLattice::PolyakovPosToIndex(x + 2, y, z)]));
        avgcound++;
      }
    }
  }
  Result /= avgcound;
  printf("%f	", Result);
  printf("\n");
  delete[] PolyakovLoops;
  return;
}

//Run Polyakov step of MC simulation
void SUNLattice::PolyakovRun() {
  int thermN = 1000;
  int N = 100;
  for (int i = 0; i < thermN; i++) {
    SUNLattice::MonteCarloUpdate();
    printf("Thermalising : %i / %i\n", i + 1, thermN);
  }
  for (int i = 0; i < N; i++) {
    SUNLattice::MonteCarloUpdate();
    SUNLattice::PolyakovLoopCorrelation();
  }
}
//Run Wilsonloop step of MC simulation
void SUNLattice::WilsonLoopRun(const char *filename, int autocorrTime, int N,
                               const vector<pair<int, int>> Loops) {

  FILE *f = fopen(filename, "a");
  fprintf(f, "## Beta : %f	L : %i	T : %i\n", settings.beta,
          settings.s_LattSize, settings.t_LattSize);

  fprintf(f, "i");
  for (int i = 0; i < static_cast<int>(Loops.size()); i++) {
    fprintf(f, "	%ix%i", Loops[i].first, Loops[i].second);
  }
  fprintf(f, "\n");
  fclose(f);
  Thermalize(20000, 0);

  for (int i = 0; i < N; i++) {
    Thermalize(autocorrTime, 0);
    f = fopen(filename, "a");
    fprintf(f, "\n%i", i);
    for (int i = 0; i < static_cast<int>(Loops.size()); i++) {
      fprintf(f, "	%f", WilsonLoop(Loops[i].first, Loops[i].second));
    }
    fclose(f);
  }
}

//Run full simulation of plaque-plaque correlation simulation
void SUNLattice::PlaqPlaqRun() {
  int thermN = 2000;
  int N = 2000000;
  for (int i = 0; i < thermN; i++) {
    SUNLattice::MonteCarloUpdate();
    fprintf(stderr, "#Thermalising : %i / %i\n", i + 1, thermN);
  }
  Complex res;
  for (int i = 0; i < N; i++) {
    printf("%i", i);
    SUNLattice::MonteCarloUpdate();
    for (int sep = 1; sep < 15; sep++) {
      res = PlaqPlaqCorrelation(sep);
      printf("	%f	%f", real(res), imag(res));
    }

    printf("\n");
  }
}

//Bring lattice into thermal equilibrium
void SUNLattice::Thermalize(int Ntherm, int verbose) {
  for (int i = 0; i < Ntherm; i++) {
    SUNLattice::MonteCarloUpdate();
    if (verbose) {
      fprintf(stdout, "\033[2J\033[0;1f");
      printf("Thermalising : %i / %i\n", i + 1, Ntherm);
    }
  }
}

//Calculate Autocorrelationtime of Wilsons
int SUNLattice::AutoCorrelationWilson(int N, int a, int b) {
  vector<double> rho;
  vector<vector<double>> autocorrInitialState(settings.Volume() * 4,
                                              vector<double>(4));
  double autocorrInitialAVG;
  vector<double> autocorrFinalState(settings.Volume() * 4);
  double autocorrFinalAVG;
  double C0, C1;
  C0 = 0.0;
  forlatt {
    LatticePosition pos(t, x, y, z, settings);
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = 0; nu < 4; nu++) {
        if (nu == mu)
          continue;
        autocorrInitialState[LatticePositionToIndex(pos, mu)][nu] =
            real(Plaquette(pos, mu, nu, a, b).trace()) / 3;
        C0 += autocorrInitialState[LatticePositionToIndex(pos, mu)][nu] *
              autocorrInitialState[LatticePositionToIndex(pos, mu)][nu];
      }
    }
  }
  C0 /= settings.Volume() * 12;
  autocorrInitialAVG = WilsonLoop(a, b);
  C0 -= autocorrInitialAVG * autocorrInitialAVG;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < 1; j++) {
      SUNLattice::MonteCarloUpdate();
    }

    C1 = 0.0;
    forlatt {
      LatticePosition pos(t, x, y, z, settings);
      for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
          if (nu == mu)
            continue;
          C1 += autocorrInitialState[LatticePositionToIndex(pos, mu)][nu] *
                real(Plaquette(pos, mu, nu, a, b).trace()) / 3;
        }
      }
    }
    C1 /= settings.Volume() * 12;
    autocorrFinalAVG = WilsonLoop(a, b);
    C1 -= autocorrInitialAVG * autocorrFinalAVG;
    rho.push_back(C1 / C0);
  }
  double t_int = 0.5;
  FILE *f = fopen("beta8AutoCorr", "a");
  fprintf(f, "## Beta : %f	L : %i	T : %i\n", settings.beta,
          settings.s_LattSize, settings.t_LattSize);
  for (int i = 0; i < N; i++) {
    fprintf(f, "%i	%f\n", i, rho[i]);
    t_int += rho[i];
  }
  fclose(f);

  return (int)t_int + 1;
}

//Run full simulation of autocorrelations of Wilsonloops (1,1)
void SUNLattice::Run(const char *filename, int autocorrTime) {
  Thermalize(20000, 1);
  printf("%2.2f	%i\n", settings.beta, AutoCorrelationWilson(1000, 1, 1));
  return;
}
}
