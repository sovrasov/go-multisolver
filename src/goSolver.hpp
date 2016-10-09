#ifndef GO_SOLVER_HPP
#define GO_SOLVER_HPP

#include "problemsPool.hpp"

struct SolverParameters
{
  double eps;
  double r;
  unsigned numThreads;
  unsigned iterationsLimit;

  SolverParameters() {}
  SolverParameters(double _eps, double _r,
      unsigned _numThreads, unsigned _iterationsLimit) :
        eps(_eps), r(_r), numThreads(_numThreads), iterationsLimit(_iterationsLimit)
  {}
};

template <class FType>
class GOSolver
{
protected:

  ProblemsPool<FType> mProblems;

public:

  void SetProblemsPool(ProblemsPool<FType>& problems);
  void SetParameters(SolverParameters& params);
  void Solve();
  GetOptimumEstimations() {}

};

template <class FType>
void GOSolver<FType>::SetProblemsPool(ProblemsPool<FType>& problems)
{

}

template <class FType>
void GOSolver<FType>::SetParameters(SolverParameters& params)
{

}

template <class FType>
void GOSolver<FType>::Solve()
{

}

#endif
