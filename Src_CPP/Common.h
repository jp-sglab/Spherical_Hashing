#pragma once

#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <bitset>

#include "Parameters.h"
#include "Stopwatch.hpp"

using namespace std;

#define		PI				3.1415926536

int Rand_NDigits(int nDigits);
int Rand_Uniform_Int(int minVal, int maxVal);

template<typename RealT>
RealT Rand_Uniform(RealT minVal, RealT maxVal)
{
	return ( minVal + ( maxVal - minVal ) * ( (RealT)rand() / (RealT)(RAND_MAX) ) );
}

template<typename RealT>
RealT Rand_Gaussian()
{
	RealT x1, x2, ret;
	do
	{
		x1 = Rand_Uniform<RealT>( 0.0 , 1.0 );
	} while( x1 == 0 );
	x2 = Rand_Uniform<RealT>( 0.0 , 1.0 );
	ret = sqrt( -2.0 * log( x1 ) ) * cos( 2.0 * PI * x2 );
	return ret;
}

template<typename RealT>
RealT Compute_Distance_L2Sq(RealT *v0, RealT *v1, int dim)
{
	RealT ret = 0.0;
	for(int i=0;i<dim;i++)
	{
		ret += ( v0[i] - v1[i] ) * ( v0[i] - v1[i] );
	}
	return ret;
}

template<typename RealT>
void SetVector_Val(RealT *vec, int dim, RealT val)
{
	for(int i=0;i<dim;i++)
	{
		vec[i] = val;
	}
}

template<typename RealT>
void SetVector_Vec(RealT *dest, RealT *src, int dim)
{
	for(int i=0;i<dim;i++)
	{
		dest[i] = src[i];
	}
}

template<typename RealT>
void Scalar_Vector(RealT *A, RealT s, int dim)
{
	for(int i=0;i<dim;i++)
	{
		A[i] = s * A[i];
	}
}

template<typename RealT>
void Sub_Vector(RealT *A, RealT *B, RealT *ret, int dim)
{
	for(int i=0;i<dim;i++)
	{
		ret[i] = A[i] - B[i];
	}
}


template<typename RealT>
void Add_Vector(RealT *A, RealT *B, RealT *ret, int dim)
{
	for(int i=0;i<dim;i++)
	{
		ret[i] = A[i] + B[i];
	}
}

template<typename RealT>
void SetMatrix_Val(RealT **mat, int nY, int nX, RealT val)
{
	for(int i=0;i<nY;i++)
	{
		for(int j=0;j<nX;j++)
		{
			mat[i][j] = val;
		}
	}
}