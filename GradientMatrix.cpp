//
//  GradientMatrix.cpp
//
//
//

#include "GradientMatrix.hpp"

GradientMatrix::GradientMatrix()
{
    //create();
}

GradientMatrix::GradientMatrix(const double a_size, const double b_size, const double c_size, const int n_num) : a(a_size), b(b_size), c(c_size), n(n_num)
{
    create();
}

GradientMatrix::~GradientMatrix()
{
    //destroy();
}

void GradientMatrix::destroy()
{
    MatDestroy(&matrix); 
    MatDestroy(&transpose);
    return 0;
}

void GradientMatrix::create()
{
    double temp;

    double tempA = a, tempB = b, tempC = c;
    
    // now the gradient (or B) matrices
    MatCreateSeqDense(PETSC_COMM_SELF, 6, 24, PETSC_NULL, &matrix);
	
	double x[8],y[8],z[8];
    double a, b, c, d, e, f, g, h, l;
    double difxs, difys, difzs, difxt, difyt, difzt, difxn, difyn, difzn;
    for (int i = 0; i < 8; i++)
    {
        x[i] = (1 + SIGN(integMtx(i, 0))) / 2.0 * tempA;
        y[i] = (1 + SIGN(integMtx(i, 1))) / 2.0 * tempB;
        z[i] = (1 + SIGN(integMtx(i, 2))) / 2.0 * tempC;
    }
    difxs = difxt = difxn = difys = difyt = difyn = difzs = difzt = difzn = 0;
    for (int i = 0; i < 8; i++)
    {
        difxs += diffMtx(n, 0, i) * x[i];
        difxt += diffMtx(n, 1, i) * x[i];
        difxn += diffMtx(n, 2, i) * x[i];
        difys += diffMtx(n, 0, i) * y[i];
        difyt += diffMtx(n, 1, i) * y[i];
        difyn += diffMtx(n, 2, i) * y[i];
        difzs += diffMtx(n, 0, i) * z[i];
        difzt += diffMtx(n, 1, i) * z[i];
        difzn += diffMtx(n, 2, i) * z[i];
    }
	
    Jdet = difxs * difyt * difzn + difxt * difyn * difzs + difxn * difys * difzt - difxn * difyt * difzs - difxt * difys * difzn - difxs * difyn * difzt;
	
    a = difyt * difzn - difzt * difyn;
    b = difys * difzn - difzs * difyn;
    c = difys * difzt - difzs * difyt;
    d = difxt * difzn - difzt * difxn;
    e = difxs * difzn - difzs * difxn;
    f = difxs * difzt - difzs * difxt;
    g = difxt * difyn - difyt * difxn;
    h = difxs * difyn - difys * difxn;
    l = difxs * difyt - difys * difxt;
	
    for (int i = 0; i < NODES_PER_ELEMENT; i++)
    {// 8 integration points

        /*tempA = 1.0 / tempA * SIGN(integMtx(i, 0)) *
            (0.5 + integMtx(n, 1) * SIGN(integMtx(i, 1))) *
            (0.5 + integMtx(n, 2) * SIGN(integMtx(i, 2)));

        tempB = 1.0 / tempB * SIGN(integMtx(i, 1)) *
            (0.5 + integMtx(n, 0) * SIGN(integMtx(i, 0))) *
            (0.5 + integMtx(n, 2) * SIGN(integMtx(i, 2)));

        tempC = 1.0 / tempC * SIGN(integMtx(i, 2)) *
            (0.5 + integMtx(n, 0) * SIGN(integMtx(i, 0))) *
            (0.5 + integMtx(n, 1) * SIGN(integMtx(i, 1)));*/

    	
        temp = a * diffMtx(n, 0, i) - b * diffMtx(n, 1, i) + c * diffMtx(n, 2, i);
        temp = temp / Jdet;
        
        MatSetValue(matrix, 0, 0 + 3 * i, temp, INSERT_VALUES);
        MatSetValue(matrix, 3, 1 + 3 * i, temp, INSERT_VALUES);
        MatSetValue(matrix, 5, 2 + 3 * i, temp, INSERT_VALUES);
        
        temp = -d * diffMtx(n, 0, i) + e * diffMtx(n, 1, i) - f * diffMtx(n, 2, i);
        temp = temp / Jdet;
        
        
        MatSetValue(matrix, 1, 1 + 3 * i, temp, INSERT_VALUES);
        MatSetValue(matrix, 3, 0 + 3 * i, temp, INSERT_VALUES);
        MatSetValue(matrix, 4, 2 + 3 * i, temp, INSERT_VALUES);
        
        temp = g * diffMtx(n, 0, i) - h * diffMtx(n, 1, i) + l * diffMtx(n, 2, i);
        temp = temp / Jdet;
        
        
        MatSetValue(matrix, 2, 2 + 3 * i, temp, INSERT_VALUES);
        MatSetValue(matrix, 4, 1 + 3 * i, temp, INSERT_VALUES);
        MatSetValue(matrix, 5, 0 + 3 * i, temp, INSERT_VALUES);
    }
	
    MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);

    // create transpose
    MatTranspose(matrix, MAT_INITIAL_MATRIX, &transpose);

    return 0;
}