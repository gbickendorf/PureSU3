#include "TwoLinkOperator.h"

TwoLinkOperator::TwoLinkOperator(matrix3 U1, matrix3 U2)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			u1.M[i][j]= U1.M[i][j];
			u2.M[i][j]= U2.M[i][j];
		}
	}
}
TwoLinkOperator::TwoLinkOperator()
{
	TwoLinkOperator(matrix3(),matrix3());
}

void TwoLinkOperator::DivideByScalar(double s)
{
		for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			u1.M[i][j]/=s;
			u2.M[i][j]/=s;
		}
	}
}

Complex TwoLinkOperator::trace()
{
	return conj(u1.trace())*u2.trace();
}

TwoLinkOperator operator*(TwoLinkOperator t1, TwoLinkOperator t2)
{
	return TwoLinkOperator(t1.u1*t2.u1,t1.u2*t2.u2);
	
}

TwoLinkOperator operator+(TwoLinkOperator t1, TwoLinkOperator t2)
{
	return TwoLinkOperator(t1.u1+t2.u1,t1.u2+t2.u2);
	
}

