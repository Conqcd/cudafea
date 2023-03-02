//
//  GradientMatrix.h 
//

#ifndef GRADIENTMATRIX_H
#define GRADIENTMATRIX_H

#include "Common.h"
#include "IntegrationMatrix.h"


class GradientMatrix
{
    
public:
    GradientMatrix();
    GradientMatrix(double a_size, double b_size, double c_size, int n_num);
    ~GradientMatrix();
    
    // get matrix and transpose
    Mat* getmat(){ return (&matrix); }
    Mat* getmat_t(){ return (&transpose); }
    double getJdet() const { return Jdet; }
	

protected:
    
    int n; // n
    double a, b, c; // voxel sizes (x, y, z)
    
    IntegrationMatrix integMtx;
    Diff_Ni_j diffMtx;
    double Jdet;
    Mat matrix; // petsc matrix
    Mat transpose; // transpose of pmatrix
    
    
    void create();
    void destroy();

};

#endif /* GRADIENTMATRIX_H */
