//
//  LocalStiffnessMatrix.h
//

#ifndef LOCALSTIFFNESSMATRIX_H
#define LOCALSTIFFNESSMATRIX_H

#include "Common.h"
#include "GradientMatrix.h"
#include "PropertyMatrix.h"


class LocalStiffnessMatrix
{
public:
    double youngsm; // youngs modulus
    double poissonsr; // poissons ratio
    double a, b, c; // voxel sizes (x, y, z)
    
    PetscErrorCode ierr; // petsc errorcode
    
    PetscScalar *matrix; // local stiffness matrix
    
    LocalStiffnessMatrix();
    LocalStiffnessMatrix(const double a, const double b, const double c, const double ym, const double pr);
    LocalStiffnessMatrix(const double a, const double b, const double c, const double ym, const double pr, GradientMatrix * GradientMatrices);
    ~LocalStiffnessMatrix();
    
    // inline
    // get matrix
    PetscScalar* getmat(){ return (matrix); }
    
    double getval(const unsigned int index)
    {
        return matrix[index];
    }
    
    
    
protected:
    
    PetscErrorCode create(GradientMatrix *GradientMatrices);
    PetscErrorCode create();
    void destroy();

};

#endif /* LOCALSTIFFNESSMATRIX_H */
