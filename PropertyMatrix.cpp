//
//  PropertyMatrix.cpp
//

#include "PropertyMatrix.hpp"

// use default values
PropertyMatrix::PropertyMatrix() : youngsm(1), poissonsr(0.3)
{
    create();
}

PropertyMatrix::PropertyMatrix(const double ym, const double pr) : youngsm(ym), poissonsr(pr)
{
    create();
}

PropertyMatrix::~PropertyMatrix()
{
    destroy();
}

void PropertyMatrix::destroy()
{
    // MatDestroy(&matrix);
}

// Build Property Matrix
void PropertyMatrix::create()
{
    matrix.reset(6,6);
    
    // MatSetValue(matrix, 0, 0, 1.0, INSERT_VALUES);
    // MatSetValue(matrix, 0, 1, poissonsr / (1.0 - poissonsr), INSERT_VALUES);
    // MatSetValue(matrix, 0, 2, poissonsr / (1.0 - poissonsr), INSERT_VALUES);
    // MatSetValue(matrix, 1, 0, poissonsr / (1.0 - poissonsr), INSERT_VALUES);
    // MatSetValue(matrix, 1, 1, 1.0, INSERT_VALUES);
    // MatSetValue(matrix, 1, 2, poissonsr / (1.0 - poissonsr), INSERT_VALUES);
    // MatSetValue(matrix, 2, 0, poissonsr / (1.0 - poissonsr), INSERT_VALUES);
    // MatSetValue(matrix, 2, 1, poissonsr / (1.0 - poissonsr), INSERT_VALUES);
    // MatSetValue(matrix, 2, 2, 1.0, INSERT_VALUES);
    // MatSetValue(matrix, 3, 3, (1.0 - 2.0 * poissonsr) / (2.0 * (1.0 - poissonsr)), INSERT_VALUES);
    // MatSetValue(matrix, 4, 4, (1.0 - 2.0 * poissonsr) / (2.0 * (1.0 - poissonsr)), INSERT_VALUES);
    // MatSetValue(matrix, 5, 5, (1.0 - 2.0 * poissonsr) / (2.0 * (1.0 - poissonsr)), INSERT_VALUES);
    
    matrix.insert(0,0,1);
    
    matrix.scale(youngsm * (1.0 - poissonsr)/((1.0 + poissonsr) * (1.0 - 2.0 * poissonsr)));
}