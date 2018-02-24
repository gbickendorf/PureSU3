#ifndef SETTING2_H
#define SETTING2_H

using namespace std;


enum class StartCondition {Hot, Cold};
 
class LatticeSettings
{
public:
	LatticeSettings(double Beta, int S_LattSize, int T_LattSize, StartCondition condition);
	double beta;
	int s_LattSize;
	int t_LattSize;
	int Volume(void);
	StartCondition Condition;
};

#endif
