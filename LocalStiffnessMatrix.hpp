//
//  LocalStiffnessMatrix.h
//

#ifndef LOCALSTIFFNESSMATRIX_H
#define LOCALSTIFFNESSMATRIX_H

#include "Common.hpp"
#include "GradientMatrix.hpp"
#include "PropertyMatrix.hpp"


class LocalStiffnessMatrix
{
public:
    double youngsm; // youngs modulus
    double poissonsr; // poissons ratio
    double a, b, c; // voxel sizes (x, y, z)
    
    
    DenseMatrix matrix; // local stiffness matrix
    
    LocalStiffnessMatrix();
    LocalStiffnessMatrix(const double a, const double b, const double c, const double ym, const double pr);
    LocalStiffnessMatrix(const double a, const double b, const double c, const double ym, const double pr, const std::vector<GradientMatrix>& GradientMatrices);
    ~LocalStiffnessMatrix();
    
    // inline
    // get matrix
    const Matrix& getmat() const{ return matrix; }
    
    double getval(const unsigned int index) const
    {
        return matrix[index];
    }
    
protected:
    
    void create(const std::vector<GradientMatrix>& GradientMatrices);
    void create();
    void destroy();
};

#endif /* LOCALSTIFFNESSMATRIX_H */
