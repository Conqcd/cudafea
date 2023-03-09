#include "BiCGSTABSolver.hpp"

extern void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b);

BiCGSTABSolver::BiCGSTABSolver()
        :iteration_set(1000),iteration_close(0),ResdualNorm(0),tolerance(1e-9)
{
}

BiCGSTABSolver::~BiCGSTABSolver()
{
}

void BiCGSTABSolver::Solve(const SymetrixSparseMatrix& A,Vector& x,const Vector& b)
{
    BiCGSTAB(A,x,b);
}

void BiCGSTABSolver::setMaxIteration(idxType iter)
{
    iteration_set = iter;
}

void BiCGSTABSolver::setTolerance(double tol)
{
    tolerance = tol;
}

