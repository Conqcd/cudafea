//
//  IntegrationMatrix.h
//  
//
//

#ifndef INTEGRATIONMATRIX_H
#define INTEGRATIONMATRIX_H

#include <vector>

#include "Common.h"

class IntegrationMatrix
{

public:
    int type;
    IntegrationMatrix() { type = 0; create(); }
    IntegrationMatrix(int _type);
    ~IntegrationMatrix();

    void create();
    
    double operator()(const unsigned int idxcol, const unsigned int idxrow) const
    { return  matrix[idxcol][idxrow]; }

private:
    double matrix[8][3];
};

class Diff_Ni_j
{

public:
    int type;
    Diff_Ni_j();
    ~Diff_Ni_j();

    void create();

    double operator()(const unsigned int idxcol, const unsigned int idxrow, const unsigned int iddepth) const
    {
        return  matrix[idxcol][idxrow][iddepth];
    }

private:

    double matrix[SAMPLES_PER_ELEMENT][DOF_3D][NODES_PER_ELEMENT];
    std::vector<double> weight;
};

#endif /* INTEGRATIONMATRIX_H */