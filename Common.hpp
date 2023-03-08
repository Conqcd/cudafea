//
//  Common.h
//  
#pragma once

 
#define SIGN(x) ((x > 0) - (x < 0))

#define DOF_3D              3   // Number of degrees of freedom per node
#define NUM_NEIGHBRS        27  // Max number of neighbours (element or nodes)
#define NUM_TERMS           81  // DOF_3D x NUM_NEIGHBRS 3 * 27
#define NODES_PER_ELEMENT   8
#define ELEMENTS_PER_NODE   8
#define SAMPLES_PER_ELEMENT 8

typedef unsigned int        xyzType;
typedef unsigned int        midxType;
// typedef uint64_t            idxType;
typedef unsigned long long  idxType;
typedef double Scalar;

#include "Math/Matrix.hpp"
#include "Math/Vector.hpp"
