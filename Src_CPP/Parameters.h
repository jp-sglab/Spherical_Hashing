#pragma once

#define	REAL_TYPE				float


// target number of nearest neighbors
#define K						1000

// binary code length
#define BCODE_LEN				64

// number of training samples for spherical hashing
#define NUM_TRAIN_SAMPLES		20000

// desired portion of training set inside of one hyper-sphere
#define INCLUDING_RATIO			0.5
// desired portion of training set inside of two hyper-spheres
#define OVERLAP_RATIO			0.25

// e_m and e_s
#define EPSILON_MEAN			0.10
#define EPSILON_STDDEV			0.15

#define MAX_NUM_ITERATIONS		50

#define INPUT_DATA_FILE_NAME	"F:\\KNN_Datasets\\BigANN_GIST_1M_960Dim\\Data.points"
#define INPUT_QUERY_FILE_NAME	"F:\\KNN_Datasets\\BigANN_GIST_1M_960Dim\\Query.points"

// to disable parallelization, comment out this
#define USE_PARALLELIZATION		