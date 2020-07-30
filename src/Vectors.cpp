#include "Vectors.h"

LatticePosition::LatticePosition(int t, int x, int y, int z,
                                 LatticeSettings &Settings)
    : settings(Settings) {
  Pos[0] = LatticePosition::wrapt(t);
  Pos[1] = LatticePosition::wraps(x);
  Pos[2] = LatticePosition::wraps(y);
  Pos[3] = LatticePosition::wraps(z);
}

//Defines periodic vector addition
int LatticePosition::wrapt(int t) {
  return (t + settings.t_LattSize) % settings.t_LattSize;
}

//Defines periodic vector addition
int LatticePosition::wraps(int s) {
  return (s + settings.s_LattSize) % settings.s_LattSize;
}

//Define periodic scalar vector multiplication
LatticePosition operator*(LatticePosition &p1, int n) {
  return LatticePosition(p1.Pos[0] * n, p1.Pos[1] * n, p1.Pos[2] * n,
                         p1.Pos[3] * n, p1.settings);
}

//Define periodic vector addition
LatticePosition operator+(LatticePosition p1, LatticePosition p2) {
  return LatticePosition(p1.Pos[0] + p2.Pos[0], p1.Pos[1] + p2.Pos[1],
                         p1.Pos[2] + p2.Pos[2], p1.Pos[3] + p2.Pos[3],
                         p1.settings);
}

//Define periodic vector subtraction 
LatticePosition operator-(LatticePosition p1, LatticePosition p2) {
  return LatticePosition(p1.Pos[0] - p2.Pos[0], p1.Pos[1] - p2.Pos[1],
                         p1.Pos[2] - p2.Pos[2], p1.Pos[3] - p2.Pos[3],
                         p1.settings);
}
