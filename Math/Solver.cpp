#include "Solver.hpp"
#include <cmath>

extern void CG(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm);
extern void PCG_ICC(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm);
extern void PCG_SSOR(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm);
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
    // CG(A,x,b,tolerance,iteration_set,iteration_close,ResdualNorm);
    // PCG_SSOR(A,x,b,tolerance,iteration_set,iteration_close,ResdualNorm);
    PCG_ICC(A,x,b,tolerance,iteration_set,iteration_close,ResdualNorm);
    
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

namespace Math
{
    
Vector solve3equation(Scalar a,Scalar b,Scalar c,Scalar d)
{
    Scalar C,B,A,Delta;
    C = c * c - 3 * b * d; 
    B = b * c - 9 * a * d;
    A = b * b - 3 * a * c;
    if(A == 0 && A == B)
    {
        Vector bb(3);
        bb[0] = bb[1] = bb[2] = -b / 3 / a;
        return bb;
    }
    Delta = B * B - 4 * A * C;
    if(Delta > 0)
        return {};
    else if(Delta == 0)
    {
        auto K = A / B;
        Vector bb(3);
        bb[0] = std::max(-b / a + K,(Scalar)-0.5 * K);
        bb[1] = bb[2] = std::min(-b / a + K,(Scalar)-0.5 * K);
        return bb;
    }else
    {
        auto T = (2 * A * b - 3 * a * B) / 2 / std::sqrt(A * A * A);       
        auto theta = std::acos(T);

        Vector bb(3);
        bb[0] = (-b - 2 * std::sqrt(A) * std::cos(theta / 3)) / 3 / a;
        bb[1] = (-b + std::sqrt(A) * (std::cos(theta / 3) + std::sqrt(3) * std::sin(theta / 3))) / 3 / a;
        bb[2] = (-b + std::sqrt(A) * (std::cos(theta / 3) - std::sqrt(3) * std::sin(theta / 3))) / 3 / a;
        if(bb[1] < bb[2])
            std::swap(bb[1],bb[2]);
        if(bb[0] < bb[1])
            std::swap(bb[0],bb[1]);
        if(bb[1] < bb[2])
            std::swap(bb[1],bb[2]);
        return bb;
    }
}

} // namespace Math