#pragma once

#include "Common.h"

class Points
{
public :
	int nP;
	int dim;
	REAL_TYPE **d;
	
	void Initialize(int _nP, int _dim)
	{
		nP = _nP;		dim = _dim;
		d = new REAL_TYPE * [ nP ];
		for(int i=0;i<nP;i++)
		{
			d[i] = new REAL_TYPE [ dim ];
		}
	}

	// this function is for read point set from file
	// format:
	// (num. points) (dimensionality)
	// v0 (floats, number of elements is equal to dimensionality)
	// v1
	// ...
	void Initialize_From_File(char *fileName)
	{
		FILE *input = fopen( fileName , "rb" );
		fread( &nP , sizeof(int) , 1 , input );
		fread( &dim , sizeof(int) , 1 , input );
		d = new REAL_TYPE * [ nP ];
		for(int i=0;i<nP;i++)
		{
			d[i] = new REAL_TYPE [ dim ];
		}
		float *tmp = new float [ dim ];
		for(int i=0;i<nP;i++)
		{
			fread( tmp , sizeof(float) , dim , input );
			for(int k=0;k<dim;k++)
			{
				d[i][k] = (REAL_TYPE)( tmp[k] );
			}
		}
		fclose(input);
		delete [] tmp;
	}

	// computing center of points for zero centering
	void Compute_Center(REAL_TYPE *center)
	{
		double *tCenter = new double [dim];
		SetVector_Val<double>( tCenter , dim , 0.0 );
		for(int i=0;i<nP;i++)
		{
			for(int k=0;k<dim;k++)
			{
				tCenter[k] += d[i][k];
			}
		}
		for(int k=0;k<dim;k++)
		{
			tCenter[k] /= (double)(nP);
			center[k] = (REAL_TYPE)( tCenter[k] );
		}
		delete [] tCenter;
	}

	void ReleaseMem()
	{
		for(int i=0;i<nP;i++)
		{
			delete [] d[i];
		}
		delete [] d;
	}

	
};