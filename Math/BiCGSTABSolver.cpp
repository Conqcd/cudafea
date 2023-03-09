#include "BiCGSTABSolver.hpp"

extern void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm);

BiCGSTABSolver::BiCGSTABSolver()
        :iteration_set(1000),iteration_close(0),ResdualNorm(0),tolerance(1e-9)
{
}

BiCGSTABSolver::~BiCGSTABSolver()
{
}

void BiCGSTABSolver::Solve(const SymetrixSparseMatrix& A,Vector& x,const Vector& b)
{
    BiCGSTAB(A,x,b,tolerance,iteration_set,iteration_close,ResdualNorm);
}

void BiCGSTABSolver::setMaxIteration(idxType iter)
{
    iteration_set = iter;
}

void BiCGSTABSolver::setTolerance(double tol)
{
    tolerance = tol;
}

