#include "LatticeSettings.h"



LatticeSettings::LatticeSettings(double Beta, int S_LattSize, int T_LattSize, StartCondition condition)
{
	beta = Beta;
	s_LattSize = S_LattSize;
	t_LattSize = T_LattSize;
	Condition = condition;
}

int LatticeSettings::Volume()
{
	return t_LattSize*s_LattSize*s_LattSize*s_LattSize;
}
