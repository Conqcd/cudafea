#include "BiCGSTABSolver.hpp"

extern void BiCGSTAB();

BiCGSTABSolver::BiCGSTABSolver()
        :iteration_set(1000),iteration_close(0),ResdualNorm(0),tolerance(1e-9)
{
}

BiCGSTABSolver::~BiCGSTABSolver()
{
}

void BiCGSTABSolver::Solve(const Matrix& A,Vector& x,const Vector& b)
{
    BiCGSTAB();
}

void BiCGSTABSolver::setMaxIteration(idxType iter)
{
    iteration_set = iter;
}

void BiCGSTABSolver::setTolerance(double tol)
{
    tolerance = tol;
}

