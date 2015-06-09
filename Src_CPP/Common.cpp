#include "Common.h"

int Rand_NDigits(int nDigits)
{
	int ret = 0;
	int tmp = 1;
	for(int i=0;i<nDigits;i++)
	{
		ret += ( rand() % 10 ) * tmp;
		tmp *= 10;
	}
	return ret;
}

int Rand_Uniform_Int(int minVal, int maxVal)
{
	int diff = maxVal - minVal + 1;
	return ( minVal + ( Rand_NDigits(9) % diff ) );
}