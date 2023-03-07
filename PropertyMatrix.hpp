#pragma once

#include "Common.hpp"
#include "Math/Matrix.hpp"

class PropertyMatrix
{
    
public:
    
    PropertyMatrix();
    PropertyMatrix(const double ym, const double pr);
    ~PropertyMatrix();

    // get underlying matrix
    const Matrix& getmat() { return matrix; }
    
protected:
    double youngsm; // young's modulus
    double poissonsr; // poisson's ratio
    
    DenseMatrix matrix; // petsc matrix

    void create();
    void destroy();
};