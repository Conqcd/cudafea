#pragma once

#include "Matrix.hpp"

class Solver
{
private:
    int iteration_set,iteration_close;
    double ResdualNorm;
    double tolerance;
public:
    Solver(/* args */);
    ~Solver();
    void Solve(const SymetrixSparseMatrix&,Vector&,const Vector&);
    void setMaxIteration(idxType);
    void setTolerance(double);

    inline auto getIteration()const {return iteration_close;}
    inline auto getResdualNorm()const {return ResdualNorm;}
};


namespace Math
{
    
Vector solve3equation(Scalar a,Scalar b,Scalar c,Scalar d);

} // namespace Math
