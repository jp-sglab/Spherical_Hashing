#pragma once

template<typename DistType>
class Result_Element
{
public :
	int index;
	DistType dist;
	bool isCorrect;
	
	bool operator < (const Result_Element<DistType> &T) const
	{
		if( this->dist < T.dist )
		{
			return true;
		}
		return false;
	}
};


// compute average precision for each query
template<typename DistType>
double Compute_AP(int *gt, Result_Element<DistType> *res, int nP)
{
	double ret = 0.0;

#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nP;i++)
	{
		res[i].isCorrect = false;
	}

#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<K;i++)
	{
		res[ gt[i] ].isCorrect = true;
	}
	sort( &res[0] , &res[nP] );
	
	int curr_nR, curr_n, currPos, pos, nCorr;
	curr_nR = 0;		curr_n = 0;		currPos = 0;
	for(int i=0;i<nP;i++)
	{
		pos = i;
		nCorr = 0;
		// to consider equal hamming / spherical hamming distance
		while(true)
		{
			if( pos == nP || res[i].dist != res[pos].dist )
			{
				break;
			}
			if( res[pos].isCorrect )
			{
				curr_nR++;		nCorr++;
			}
			pos++;
		}
		ret += (double)( curr_nR ) / (double)(pos) * (double)nCorr;
		i = pos - 1;
	}
	ret /= (double)(K);
	return ret;
}