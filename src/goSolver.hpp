#ifndef GO_SOLVER_HPP
#define GO_SOLVER_HPP

#include "problemsPool.hpp"
#include "dataTypes.hpp"
#include "evolvent.hpp"
#include "intervalsQueue.hpp"

#include <vector>
#include <set>
#include <queue>
#include <algorithm>

enum StopType {Accuracy, OptimumVicinity, OptimalValue};

struct SolverParameters
{
  double eps;
  double r;
  unsigned numThreads;
  unsigned iterationsLimit;
  unsigned evloventTightness;

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

  SolverParameters mParameters;
  ProblemsPool<FType> mProblems;
  std::vector<bool> mActiveProblemsMask;
  unsigned mIterationsCounter;
  std::vector<Trial> mOptimumEstimations;
  std::vector<Interval> mNextIntervals;
  std::vector<Trial> mNextPoints;
  bool mNeeRefillQueue;
  Evolvent mEvolvent;
  IntervalsQueue mQueue;
  std::vector<std::set<Interval>> mSearchInformations;

  void InitDataStructures();
  void FirstIteration();
  bool UpdateHConsts();

public:

  void SetProblemsPool(ProblemsPool<FType>& problems);
  void SetParameters(SolverParameters& params);
  void Solve();
  std::vector<Trial> GetOptimumEstimations();

};


template <class FType>
void GOSolver<FType>::Solve()
{
  InitDataStructures();
  FirstIteration();
}

template <class FType>
void GOSolver<FType>::InitDataStructures()
{
  mQueue.Clear();
  mOptimumEstimations.resize(mProblems.GetNumberOfProblems());
  mActiveProblemsMask.resize(mProblems.GetNumberOfProblems());
  std::fill(mActiveProblemsMask.begin(), mActiveProblemsMask.end(), true);
  mNextIntervals.resize(mParameters.numThreads);
  mNextPoints.resize(mParameters.numThreads);
  mEvolvent = Evolvent(mProblems.GetDimension(), 12);

  Interval firstInterval;
  firstInterval.xl = 0.0;
  firstInterval.xr = 1.0;
  firstInterval.delta = 1.0;

  for (size_t i = 0; i < mProblems.GetNumberOfProblems(); i++)
  {
    double y[solverMaxDim];
    firstInterval.problemIdx = i;
    mEvolvent.GetImage(firstInterval.xl, y);
    firstInterval.zl = mProblems.CalculateObjective(y, i);
    mEvolvent.GetImage(firstInterval.xr, y);
    firstInterval.zr = mProblems.CalculateObjective(y, i);
    mSearchInformations[i].insert(firstInterval);
    mNextPoints[i].x = 0.5;
  }
}

template <class FType>
void GOSolver<FType>::FirstIteration()
{
  mIterationsCounter = 1;
}

template <class FType>
void GOSolver<FType>::SetProblemsPool(ProblemsPool<FType>& problems)
{
  mProblems = problems;
}

template <class FType>
void GOSolver<FType>::SetParameters(SolverParameters& params)
{
  mParameters = params;
}

template <class FType>
std::vector<Trial> GOSolver<FType>::GetOptimumEstimations()
{
  return mOptimumEstimations;
}

#endif
