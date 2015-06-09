#pragma once

#include "Common.h"
#include "Points.h"


// hamming distance function
__forceinline int Compute_HD(bitset<BCODE_LEN> &a, bitset<BCODE_LEN> &b)
{
	return ( ( a ^ b ).count() );
}

// spherical hamming distance function
__forceinline double Compute_SHD(bitset<BCODE_LEN> &a, bitset<BCODE_LEN> &b)
{
	return ( ( (double)( ( a ^ b ).count() ) ) / ( (double)( ( a & b ).count() ) + 0.1 ) );
}

class LSH
{
public :
	int dim;
	REAL_TYPE **pM;

	void Initialize(int _dim)
	{
		dim = _dim;
		pM = new REAL_TYPE * [ dim ];
		for(int k=0;k<dim;k++)
		{
			pM[k] = new REAL_TYPE [ BCODE_LEN ];
			for(int i=0;i<BCODE_LEN;i++)
			{
				pM[k][i] = Rand_Gaussian<REAL_TYPE>();
			}
		}
	}

	__forceinline void Compute_BCode(REAL_TYPE *x, bitset<BCODE_LEN> &y)
	{
		REAL_TYPE tmp;
		for(int i=0;i<BCODE_LEN;i++)
		{
			tmp = 0.0;
			for(int k=0;k<dim;k++)
			{
				tmp += x[k] * pM[k][i];
			}
			if( tmp > 0.0 )
			{
				y[i] = 1;
			}
			else
			{
				y[i] = 0;
			}
		}
	}
};

class Index_Distance
{
public :
	bool flag;
	int index;
	REAL_TYPE dist, distSq;
	bool operator < (const Index_Distance &T) const
	{
		if( this->distSq < T.distSq )	{			return true;		}
		return false;
	}
};

class Sphere
{
public :
	REAL_TYPE *c, r, rSq;
	
	void Initialize(int _dim)
	{
		c = new REAL_TYPE [ _dim ];
		r = 0.0;		rSq = 0.0;
	}

	// function to set radius to include desired portion of training set
	void Set_Radius(Points *ps, Index_Distance *ids);
};

class SphericalHashing
{
public :
	Points *ps;

	// training set
	Points tps;

	Sphere *s;

	Index_Distance **ids;
	bitset<NUM_TRAIN_SAMPLES> *table;

	void Initialize(Points *_ps);
	void Compute_Table();
	void Compute_Num_Overlaps(int **overlaps);
	void Set_Spheres();

	void ReleaseMem();

	__forceinline void Compute_BCode(REAL_TYPE *x, bitset<BCODE_LEN> &y)
	{
		for(int i=0;i<BCODE_LEN;i++)
		{
			if( Compute_Distance_L2Sq<REAL_TYPE>( s[i].c , x , ps->dim ) > s[i].rSq )
			{
				y[i] = 0;
			}
			else
			{
				y[i] = 1;
			}
		}
	}
};