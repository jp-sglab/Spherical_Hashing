#include "BinaryHash.h"
#include "Evaluation.h"

#ifdef USE_PARALLELIZATION
#include <omp.h>
#endif

// dps: data points set
// qps: query points set
Points dps, qps;


LSH lsh;
SphericalHashing sh;

// nP: number of data points
// nQ: number of query points
int nP, nQ;

// gt: ground truth
int **gt;

REAL_TYPE *dataCenter;


// initialize data and query points
void Initialize_Data()
{
#ifdef INPUT_DATA_FILE_NAME
	dps.Initialize_From_File( INPUT_DATA_FILE_NAME );
#endif

#ifdef INPUT_QUERY_FILE_NAME
	qps.Initialize_From_File( INPUT_QUERY_FILE_NAME );
#endif

	nP = dps.nP;		nQ = qps.nP;

	// you can control the number of queries here
	nQ = 100;

	dataCenter = new REAL_TYPE[ dps.dim ];
	// compute mean position of data points
	dps.Compute_Center( dataCenter );
}

// translate both data points and query points 
void Do_ZeroCentering()
{
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nP;i++)
	{
		Sub_Vector<REAL_TYPE>( dps.d[i] , dataCenter , dps.d[i] , dps.dim );
	}
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nQ;i++)
	{
		Sub_Vector<REAL_TYPE>( qps.d[i] , dataCenter , qps.d[i] , qps.dim );
	}
}

// undo zero-centering
void Undo_ZeroCentering()
{
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nP;i++)
	{
		Add_Vector<REAL_TYPE>( dps.d[i] , dataCenter , dps.d[i] , dps.dim );
	}
#ifdef USE_PARALLELIZATIO
	#pragma omp parallel for
#endif
	for(int i=0;i<nQ;i++)
	{
		Add_Vector<REAL_TYPE>( qps.d[i] , dataCenter , qps.d[i] , qps.dim );
	}
}

// compute ground-truth (very time comsuming)
// i recommand you to store ground truth information to the file 
void Compute_GroundTruth()
{
	gt = new int * [ nQ ];
	for(int i=0;i<nQ;i++)
	{
		gt[i] = new int [ K ];
	}
	Result_Element<REAL_TYPE> *tmp = new Result_Element<REAL_TYPE>[ nP ];
	for(int i=0;i<nQ;i++)
	{
#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
#endif
		for(int j=0;j<nP;j++)
		{
			tmp[j].index = j;
			tmp[j].dist = Compute_Distance_L2Sq<REAL_TYPE>( qps.d[i] , dps.d[j] , dps.dim );
		}
		sort( &tmp[0] , &tmp[nP] );
		#pragma omp parallel for
		for(int j=0;j<K;j++)
		{
			gt[i][j] = tmp[j].index;
		}
	}
	delete [] tmp;
}

void Process()
{
	Stopwatch T0("");

	lsh.Initialize( dps.dim );
	sh.Initialize( &dps );

	T0.Reset();		T0.Start();
	sh.Set_Spheres();
	T0.Stop();
	printf("- Learning Spherical Hashing Finished (%f seconds)\n",T0.GetTime());

	bitset<BCODE_LEN> *bCodeData_LSH = new bitset<BCODE_LEN> [ nP ];
	bitset<BCODE_LEN> *bCodeQuery_LSH = new bitset<BCODE_LEN> [ nQ ];
	bitset<BCODE_LEN> *bCodeData_SH = new bitset<BCODE_LEN> [ nP ];
	bitset<BCODE_LEN> *bCodeQuery_SH = new bitset<BCODE_LEN> [ nQ ];

	Result_Element<int> *resLSH = new Result_Element<int> [ nP ];
	Result_Element<int> *resSH_HD = new Result_Element<int> [ nP ];
	Result_Element<double> *resSH_SHD = new Result_Element<double> [ nP ];
	
	Do_ZeroCentering();

	T0.Reset();		T0.Start();
	// compute binary codes of LSH
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nP;i++)
	{
		lsh.Compute_BCode( dps.d[i] , bCodeData_LSH[i] );
		
	}
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nQ;i++)
	{
		lsh.Compute_BCode( qps.d[i] , bCodeQuery_LSH[i] );
		
	}
	T0.Stop();
	printf("- LSH: Computing Binary Codes Finished (%f seconds)\n",T0.GetTime() );

	Undo_ZeroCentering();


	T0.Reset();		T0.Start();
	// compute binary codes of Spherical Hashing
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nP;i++)
	{
		sh.Compute_BCode( dps.d[i] , bCodeData_SH[i] );
	}

#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nQ;i++)
	{
		sh.Compute_BCode( qps.d[i] , bCodeQuery_SH[i] );
	}
	T0.Stop();
	printf("- Spherical Hashing: Computing Binary Codes Finished (%f seconds)\n",T0.GetTime() );

	double mAP_LSH, mAP_SH_HD, mAP_SH_SHD;
	mAP_LSH = 0.0;		mAP_SH_HD = 0.0;		mAP_SH_SHD = 0.0;
	// process queries
	for(int qIndex=0;qIndex<nQ;qIndex++)
	{
#ifdef USE_PARALLELIZATION		
		#pragma omp parallel for
#endif
		for(int i=0;i<nP;i++)
		{
			resLSH[i].index    = i;		resLSH[i].dist    = Compute_HD(  bCodeQuery_LSH[qIndex] , bCodeData_LSH[i] );
			resSH_HD[i].index  = i;		resSH_HD[i].dist  = Compute_HD(  bCodeQuery_SH[qIndex]  , bCodeData_SH[i]  );
			resSH_SHD[i].index = i;		resSH_SHD[i].dist = Compute_SHD( bCodeQuery_SH[qIndex]  , bCodeData_SH[i]  );
		}
		mAP_LSH += Compute_AP<int>( gt[qIndex] , resLSH , nP );
		mAP_SH_HD += Compute_AP<int>( gt[qIndex] , resSH_HD , nP );
		mAP_SH_SHD += Compute_AP<double>( gt[qIndex] , resSH_SHD , nP );
	}
	mAP_LSH    /= (double)(nQ);
	mAP_SH_HD  /= (double)(nQ);
	mAP_SH_SHD /= (double)(nQ);

	printf("\n");
	printf("-- mAP\n");
	printf("\tLocality Sensitive Hashing : %f\n",mAP_LSH   );
	printf("\t    Spherical Hashing (HD) : %f\n",mAP_SH_HD );
	printf("\t    Spherical Hashing (SHD): %f\n",mAP_SH_SHD);
}

int main()
{	
	srand( (unsigned int)( time(NULL) ) );

	Stopwatch T0("");
	T0.Reset();		T0.Start();
	Initialize_Data();
	T0.Stop();
	printf("- Reading Data Finished (%f seconds)\n",T0.GetTime() );

	T0.Reset();		T0.Start();
	Compute_GroundTruth();
	T0.Stop();
	printf("- Computing GroundTruth Finished (%f seconds)\n",T0.GetTime());
	
	Process();
	
	return 1;
}