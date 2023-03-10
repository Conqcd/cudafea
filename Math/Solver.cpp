#include "Solver.hpp"

extern void CG(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm);
extern void PCG(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm);
extern void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm);

Solver::Solver()
        :iteration_set(1000),iteration_close(0),ResdualNorm(0),tolerance(1e-9)
{
}

Solver::~Solver()
{
}

void Solver::Solve(const SymetrixSparseMatrix& A,Vector& x,const Vector& b)
{
    PCG(A,x,b,tolerance,iteration_set,iteration_close,ResdualNorm);
    // BiCGSTAB(A,x,b,tolerance,iteration_set,iteration_close,ResdualNorm);
}

void Solver::setMaxIteration(idxType iter)
{
    iteration_set = iter;
}

void Solver::setTolerance(double tol)
{
    tolerance = tol;
}

