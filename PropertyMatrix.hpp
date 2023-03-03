//
//  PropertyMatrix.h
//  
//
//

#ifndef PROPERTYMATRIX_H
#define PROPERTYMATRIX_H

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
    
    Matrix matrix; // petsc matrix

    void create();
    void destroy();
};



#endif /* PROPERTYMATRIX_H */
