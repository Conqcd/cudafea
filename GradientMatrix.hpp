//
//  GradientMatrix.h 
//

#ifndef GRADIENTMATRIX_H
#define GRADIENTMATRIX_H

#include "Common.hpp"
#include "IntegrationMatrix.hpp"
#include "Math/Matrix.hpp"


class GradientMatrix
{
    
public:
    GradientMatrix();
    GradientMatrix(double a_size, double b_size, double c_size, int n_num);
    ~GradientMatrix();
    
    // get matrix and transpose
    const Matrix& getmat(){ return matrix; }
    const Matrix& getmat_t(){ return transpose; }
    double getJdet() const { return Jdet; }

protected:
    
    int n; // n
    double a, b, c; // voxel sizes (x, y, z)
    
    IntegrationMatrix integMtx;
    Diff_Ni_j diffMtx;
    double Jdet;
    Matrix matrix; // petsc matrix
    Matrix transpose; // transpose of pmatrix
    
    
    void create();
    void destroy();

};

#endif /* GRADIENTMATRIX_H */
