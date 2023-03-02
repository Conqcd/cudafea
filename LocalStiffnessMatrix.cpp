//
//  LocalStiffnessMatrix.cpp
//  
//
//

#include "LocalStiffnessMatrix.h"


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


PetscErrorCode LocalStiffnessMatrix::create()
{
    GradientMatrix gradientMtx[NODES_PER_ELEMENT];
    
    // create gradient matrices
    for(int n=0; n<NODES_PER_ELEMENT; ++n)
    {
        gradientMtx[n] = GradientMatrix(a, b, c, n);
    }
    
    // create LSM
    ierr = create(gradientMtx);
    
    return 0;
}


PetscErrorCode LocalStiffnessMatrix::create(GradientMatrix * gradientMtx)
{
        matrix = new PetscScalar();
        
        Mat LSM, propXgradMtx, tempMtx;
        int lsmlen = NODES_PER_ELEMENT * DOF_3D; // 24
        
        // Build Property Matrix
        PropertyMatrix propMtx(youngsm, poissonsr);
        
        ierr = MatCreateSeqDense(PETSC_COMM_SELF, lsmlen, lsmlen, PETSC_NULL, &LSM);
        
        // assuming all 8 gradient matrices are available
        for (int n = 0; n < SAMPLES_PER_ELEMENT; ++n)
        {// per node n
            
            ierr = MatMatMult(*(propMtx.getmat()), *(gradientMtx[n].getmat()), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &propXgradMtx); 
            ierr = MatMatMult(*(gradientMtx[n].getmat_t()), propXgradMtx, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tempMtx);
            
            ierr = MatAXPY(LSM, gradientMtx[n].getJdet(), tempMtx, DIFFERENT_NONZERO_PATTERN);
            
            ierr = MatDestroy(&propXgradMtx);
            ierr = MatDestroy(&tempMtx);
            
        }// per integration point i, gradient matrix

        //ierr = MatScale(LSM, 0.125 * a * b * c); CHKERRQ(ierr);
        //ierr = MatScale(LSM, 1.0/1.1); CHKERRQ(ierr);
        //MatView(LSM, PETSC_VIEWER_STDOUT_WORLD);
        // Get LSM
        ierr = MatDenseGetArray(LSM, (&matrix));

        return 0;
}