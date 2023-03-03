//
//  LocalStiffnessMatrix.cpp
//  
//
//

#include "LocalStiffnessMatrix.hpp"


LocalStiffnessMatrix::LocalStiffnessMatrix() : a(1.0), b(1.0), c(1.0), youngsm(1.0), poissonsr(0.3)
{}

LocalStiffnessMatrix::LocalStiffnessMatrix(const double a_size, const double b_size, const double c_size, const double ym, const double pr) : a(a_size), b(b_size), c(c_size), youngsm(ym), poissonsr(pr)
{
    create();
}

LocalStiffnessMatrix::LocalStiffnessMatrix(const double a_size, const double b_size, const double c_size, const double ym, const double pr, GradientMatrix *GradientMatrices):
a(a_size), b(b_size), c(c_size), youngsm(ym), poissonsr(pr)
{
    create(GradientMatrices);
}

LocalStiffnessMatrix::~LocalStiffnessMatrix()
{}


void LocalStiffnessMatrix::create()
{
    GradientMatrix gradientMtx[NODES_PER_ELEMENT];
    
    // create gradient matrices
    for(int n=0; n<NODES_PER_ELEMENT; ++n)
    {
        gradientMtx[n] = GradientMatrix(a, b, c, n);
    }
    
    // create LSM
    create(gradientMtx);
}


void LocalStiffnessMatrix::create(const GradientMatrix& gradientMtx)
{
    int lsmlen = NODES_PER_ELEMENT * DOF_3D; // 24
    Matrix LSM(lsmlen,lsmlen), propXgradMtx, tempMtx;
    
    // Build Property Matrix
    PropertyMatrix propMtx(youngsm, poissonsr);
    
    
    // assuming all 8 gradient matrices are available
    for (int n = 0; n < SAMPLES_PER_ELEMENT; ++n)
    {// per node n
        
        MatMatMult(*(propMtx.getmat()), *(gradientMtx[n].getmat()), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &propXgradMtx); 
        MatMatMult(*(gradientMtx[n].getmat_t()), propXgradMtx, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tempMtx);
        
        MatAXPY(LSM, gradientMtx[n].getJdet(), tempMtx, DIFFERENT_NONZERO_PATTERN);
    }// per integration point i, gradient matrix

    //ierr = MatScale(LSM, 0.125 * a * b * c); CHKERRQ(ierr);
    //ierr = MatScale(LSM, 1.0/1.1); CHKERRQ(ierr);
    //MatView(LSM, PETSC_VIEWER_STDOUT_WORLD);
    // Get LSM
    MatDenseGetArray(LSM, (&matrix));
}