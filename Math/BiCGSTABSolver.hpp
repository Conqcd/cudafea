#pragma once

#include "Matrix.hpp"

class BiCGSTABSolver
{
private:
    int iteration_set,iteration_close;
    float ResdualNorm;
    double tolerance;
public:
    BiCGSTABSolver(/* args */);
    ~BiCGSTABSolver();
    void Solve(const SymetrixSparseMatrix&,Vector&,const Vector&);
    void setMaxIteration(idxType);
    void setTolerance(double);

    inline int getIteration()const {return iteration_close;}
    inline float getResdualNorm()const {return ResdualNorm;}
};

