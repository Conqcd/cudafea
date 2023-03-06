//
//  LocalStiffnessMatrix.cpp
//  
//
//

#include "LocalStiffnessMatrix.hpp"


LocalStiffnessMatrix::LocalStiffnessMatrix() : a(1.0), b(1.0), c(1.0), youngsm(1.0), poissonsr(0.3)
{}

LocalStiffnessMatrix::LocalStiffnessMatrix(const double a_size, const double b_size, const double c_size, const double ym, const double pr) 
    :   a(a_size), b(b_size), c(c_size), youngsm(ym), poissonsr(pr)
{
    create();
}

LocalStiffnessMatrix::LocalStiffnessMatrix(const double a_size, const double b_size, const double c_size, const double ym, const double pr,const std::vector<GradientMatrix>& GradientMatrices)
    :   a(a_size), b(b_size), c(c_size), youngsm(ym), poissonsr(pr)
{
    create(GradientMatrices);
}

LocalStiffnessMatrix::~LocalStiffnessMatrix()
{}


void LocalStiffnessMatrix::create()
{
    std::vector<GradientMatrix> gradientMtx(SAMPLES_PER_ELEMENT);
    
    // create gradient matrices
    for(int n = 0; n < gradientMtx.size(); ++n)
    {
        gradientMtx[n] = GradientMatrix(a, b, c, n);
    }
    
    // create LSM
    create(gradientMtx);
}


void LocalStiffnessMatrix::create(const std::vector<GradientMatrix>& gradientMtx)
{
    int lsmlen = NODES_PER_ELEMENT * DOF_3D; // 24
    Matrix LSM(lsmlen,lsmlen), propXgradMtx, tempMtx;
    
    // Build Property Matrix
    PropertyMatrix propMtx(youngsm, poissonsr);
    
    
    // assuming all 8 gradient matrices are available
    for (int n = 0; n < gradientMtx.size(); ++n)
    {// per node n
        
        // MatMatMult(*(propMtx.getmat()), *(gradientMtx[n].getmat()), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &propXgradMtx); 
        propXgradMtx.mult(propMtx.getmat(),gradientMtx[n].getmat());
        // MatMatMult(*(gradientMtx[n].getmat_t()), propXgradMtx, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tempMtx);
        tempMtx.mult(gradientMtx[n].getmat_t(),propXgradMtx);
        
        // MatAXPY(LSM, gradientMtx[n].getJdet(), tempMtx, DIFFERENT_NONZERO_PATTERN);
        LSM.AXPY(gradientMtx[n].getJdet(),tempMtx);
    }// per integration point i, gradient matrix

    //ierr = MatScale(LSM, 0.125 * a * b * c); CHKERRQ(ierr);
    //ierr = MatScale(LSM, 1.0/1.1); CHKERRQ(ierr);
    //MatView(LSM, PETSC_VIEWER_STDOUT_WORLD);
    // Get LSM
    matrix = LSM;
    // MatDenseGetArray(LSM, (&matrix));
}