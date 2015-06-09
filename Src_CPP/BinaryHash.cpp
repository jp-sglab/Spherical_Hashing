#include "BinaryHash.h"

#include <omp.h>

void Sphere::Set_Radius(Points *ps, Index_Distance *ids)
{
	for(int i=0;i<ps->nP;i++)
	{
		ids[i].index = i;
		ids[i].distSq = Compute_Distance_L2Sq<REAL_TYPE>( c , ps->d[i] , ps->dim );
		ids[i].dist = sqrt( ids[i].distSq );
	}
	sort( &ids[0] , &ids[ ps->nP ] );
	int tIndex = (int)( (REAL_TYPE)( ps->nP ) * INCLUDING_RATIO );
	r = ids[tIndex-1].dist;
	rSq = r * r;
}

void SphericalHashing::Initialize(Points *_ps)
{
	ps = _ps;
	
	tps.Initialize( NUM_TRAIN_SAMPLES , ps->dim );	
	bool *checkList = new bool [ ps->nP ];
	SetVector_Val<bool>( checkList , ps->nP , true );
	int index;

	// sampling training set
	for(int i=0;i<NUM_TRAIN_SAMPLES;i++)
	{
		while(true)
		{
			index = Rand_Uniform_Int( 0 , ps->nP - 1 );
			if( checkList[index] )
			{
				checkList[index] = false;
				break;
			}
		}
		SetVector_Vec<REAL_TYPE>( tps.d[i] , ps->d[index] , tps.dim );		
	}
	delete [] checkList;

	s = new Sphere [ BCODE_LEN ];
	for(int i=0;i<BCODE_LEN;i++)
	{
		s[i].Initialize( tps.dim );
	}

	table = new bitset<NUM_TRAIN_SAMPLES> [ BCODE_LEN ];
	for(int i=0;i<BCODE_LEN;i++)
	{
		for(int j=0;j<NUM_TRAIN_SAMPLES;j++)
		{
			table[i][j] = 0;
		}
	}
	ids = new Index_Distance * [ BCODE_LEN ];
	for(int i=0;i<BCODE_LEN;i++)
	{
		ids[i] = new Index_Distance [ NUM_TRAIN_SAMPLES ];
	}
}

void SphericalHashing::ReleaseMem()
{
	tps.ReleaseMem();
	delete [] table;
	for(int i=0;i<BCODE_LEN;i++)
	{
		delete [] ids[i];
	}
	delete ids;
}

void SphericalHashing::Compute_Table()
{
	for(int i=0;i<BCODE_LEN;i++)
	{
#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
#endif
		for(int j=0;j<NUM_TRAIN_SAMPLES;j++)
		{
			table[i][j] = 0;
			if( Compute_Distance_L2Sq<REAL_TYPE>( s[i].c , tps.d[j] , tps.dim ) < s[i].rSq )
			{
				table[i][j] = 1;
			}
		}
	}
}

void SphericalHashing::Compute_Num_Overlaps(int **overlaps)
{
	for(int i=0;i<BCODE_LEN-1;i++)
	{
		overlaps[i][i] = table[i].count();
#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
#endif
		for(int j=i+1;j<BCODE_LEN;j++)
		{
			overlaps[i][j] = ( table[i] & table[j] ).count();
			overlaps[j][i] = overlaps[i][j];
		}
	}
}



void SphericalHashing::Set_Spheres()
{
	REAL_TYPE	marginT = 0.05;

	double allowedErrorMean, allowedErrorVar;
	allowedErrorMean = (double)tps.nP * OVERLAP_RATIO * EPSILON_MEAN;
	allowedErrorVar = (double)tps.nP * OVERLAP_RATIO * EPSILON_STDDEV;

	// initial pivots are determined as center of 10 randomly sampled points
	// for locating initial pivots near the data center
	int rIndex;
	for(int i=0;i<BCODE_LEN;i++)
	{
		SetVector_Val<REAL_TYPE>( s[i].c , tps.dim , 0.0 );
		for(int j=0;j<10;j++)
		{
			rIndex = Rand_Uniform_Int( 0 , tps.nP - 1 );
			Add_Vector<REAL_TYPE>( s[i].c , tps.d[ rIndex ] , s[i].c , tps.dim );
		}
		Scalar_Vector<REAL_TYPE>( s[i].c , 0.1 , tps.dim );
	}

#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<BCODE_LEN;i++)
	{
		s[i].Set_Radius( &tps , &ids[i][0] );
	}
	
	int **overlaps = new int * [ BCODE_LEN ];
	for(int i=0;i<BCODE_LEN;i++)
	{
		overlaps[i] = new int [ BCODE_LEN ];
	}
	REAL_TYPE **forces = new REAL_TYPE * [ BCODE_LEN ];
	for(int i=0;i<BCODE_LEN;i++)
	{
		forces[i] = new REAL_TYPE [ tps.dim ];
	}
	REAL_TYPE *force = new REAL_TYPE [ tps.dim ];
	REAL_TYPE tmpOverlap, alpha;

	for(int k=0;k<MAX_NUM_ITERATIONS;k++)
	{
		Compute_Table();
		Compute_Num_Overlaps( overlaps );

		double mean, variance, cnt;
		mean = 0.0;		cnt = 0.0;		variance = 0.0;
		for(int i=0;i<BCODE_LEN-1;i++)
		{
			for(int j=i+1;j<BCODE_LEN;j++)
			{
				mean += (double)( overlaps[i][j] );
				cnt += 1.0;
			}
		}
		mean /= cnt;
		for(int i=0;i<BCODE_LEN-1;i++)
		{
			for(int j=i+1;j<BCODE_LEN;j++)
			{
				variance += ( (double)( overlaps[i][j] ) - mean ) * ( (double)( overlaps[i][j] ) - mean );
			}
		}
		variance /= cnt;
		
		// iteration convergence condition
		if( fabs( mean - ( (double)tps.nP * OVERLAP_RATIO ) ) < allowedErrorMean && sqrt(variance) < allowedErrorVar )
		{
			break;
		}
		

		// force computation
		SetMatrix_Val<REAL_TYPE>( forces , BCODE_LEN , tps.dim , 0.0 );
		for(int i=0;i<BCODE_LEN-1;i++)
		{
			for(int j=i+1;j<BCODE_LEN;j++)
			{
				tmpOverlap = (REAL_TYPE)overlaps[i][j] / (REAL_TYPE)tps.nP;
				alpha = ( tmpOverlap - OVERLAP_RATIO ) / OVERLAP_RATIO;
				alpha /= 2.0;

				Sub_Vector<REAL_TYPE>( s[j].c , s[i].c , force , tps.dim );
				Scalar_Vector<REAL_TYPE>( force , alpha , tps.dim);
				Add_Vector<REAL_TYPE>( forces[j] , force , forces[j] , tps.dim );
				Scalar_Vector<REAL_TYPE>( force , -1.0 , tps.dim );
				Add_Vector<REAL_TYPE>( forces[i] , force , forces[i] , tps.dim );
			}
		}

		// move pivots and adjust radius
#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
#endif
		for(int i=0;i<BCODE_LEN;i++)
		{
			Scalar_Vector<REAL_TYPE>( forces[i] , 1.0 / BCODE_LEN , tps.dim );
			Add_Vector<REAL_TYPE>( s[i].c , forces[i] , s[i].c , tps.dim );
			s[i].Set_Radius( &tps , &ids[i][0] );
		}
	}

	for(int i=0;i<BCODE_LEN;i++)
	{
		delete [] overlaps[i];
		delete [] forces[i];
	}
	delete [] overlaps;
	delete [] forces;
	delete [] force;
}