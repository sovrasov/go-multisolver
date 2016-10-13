#ifndef GO_SOLVER_HPP
#define GO_SOLVER_HPP

#include "problemsPool.hpp"
#include "dataTypes.hpp"
#include "evolvent.hpp"
#include "intervalsQueue.hpp"

#include <vector>
#include <set>
#include <queue>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace solver_consts
{
  const double zeroHLevel = 1e-12;
}


enum StopType {Accuracy, OptimumVicinity, OptimalValue};

struct SolverParameters
{
  double eps;
  double r;
  unsigned numThreads;
  unsigned iterationsLimit;
  unsigned evloventTightness = 12;

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
  std::vector<double> mHEstimations;
  unsigned mIterationsCounter;
  std::vector<Trial> mOptimumEstimations;
  std::vector<Interval*> mNextIntervals;
  std::vector<Trial> mNextPoints;
  bool mNeeRefillQueue;
  Evolvent mEvolvent;
  std::priority_queue<Interval*, std::vector<Interval*>, CompareByR> mQueue;
  std::vector<std::set<Interval*>> mSearchInformations;
  unsigned mNumberOfActiveProblems;
  unsigned mNumberOfTrials;

  void InitDataStructures();
  void FirstIteration();
  void ClearDataStructures();
  void MakeTrials();
  double CalculateR(const Interval*);
  void InsertIntervals();
  void UpdateH(const Interval*);
  void EstimateOptimums();
  void RefillQueue();
  void CalculateNextPoints();
  bool CheckStopCondition();

public:

  void SetProblemsPool(ProblemsPool<FType>& problems);
  void SetParameters(SolverParameters& params);
  void Solve();
  std::vector<Trial> GetOptimumEstimations();
  unsigned GetTrialsNumber() const { return mNumberOfTrials; }
  unsigned GetIterationsNumber() const { return mIterationsCounter; }
};

template <class FType>
void GOSolver<FType>::Solve()
{
  bool needStop = false;
  InitDataStructures();
  FirstIteration();

  do {
    MakeTrials();
    InsertIntervals();//modification of refill flag
    EstimateOptimums();//modification of refill flag
    if (mNeeRefillQueue || mQueue.size() < mParameters.numThreads)
      RefillQueue();
    CalculateNextPoints();
    needStop = CheckStopCondition();
    mIterationsCounter++;
  } while(mIterationsCounter < mParameters.iterationsLimit && !needStop);

  ClearDataStructures();
}

template <class FType>
bool GOSolver<FType>::CheckStopCondition()
{
  bool needStop = false;
  for (size_t i = 0; i < mParameters.numThreads; i++)
  {
    if (mNextIntervals[i]->delta < mParameters.eps)
    {
      mActiveProblemsMask[mNextIntervals[i]->problemIdx] = false;
      mNumberOfActiveProblems--;
    }
  }

  if (mNumberOfActiveProblems == 0)
    needStop = true;

  return needStop;
}

template <class FType>
void GOSolver<FType>::CalculateNextPoints()
{
  for(size_t i = 0; i < mParameters.numThreads; i++)
  {
    mNextIntervals[i] = mQueue.top();
    mQueue.pop();
    mNextPoints[i].x = 0.5 * (mNextIntervals[i]->xl + mNextIntervals[i]->xr) -
        (((mNextIntervals[i]->zr - mNextIntervals[i]->zl) > 0) ? 1: -1) *
        pow(fabs(mNextIntervals[i]->zr - mNextIntervals[i]->zl) /
        mHEstimations[i], mProblems.GetDimension()) / 2. / mParameters.r;

    mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y);
  }
}

template <class FType>
void GOSolver<FType>::RefillQueue()
{
  mQueue = {};

  for(size_t i = 0; i < mProblems.GetNumberOfProblems(); i++)
    if(mActiveProblemsMask[i])
    {
      for(auto it = mSearchInformations[i].begin(); it != mSearchInformations[i].end(); ++it)
      {
        auto interval = *it;
        interval->R = CalculateR(interval);
        mQueue.push(interval);
      }
    }

  mNeeRefillQueue = false;
}

template <class FType>
void GOSolver<FType>::EstimateOptimums()
{
  for(size_t i = 0; i < mParameters.numThreads; i++)
  {
    unsigned problemIdx = mNextIntervals[i]->problemIdx;
    if(mNextPoints[i].z < mOptimumEstimations[problemIdx].z)
    {
      mOptimumEstimations[problemIdx] = mNextPoints[i];
      mNeeRefillQueue = true;
    }
  }
}

template <class FType>
void GOSolver<FType>::UpdateH(const Interval* i)
{
  double intervalH = fabs(i->zr - i->zl) / pow(i->xr - i->xl,
    1. / mProblems.GetDimension());
  double oldH = mHEstimations[i->problemIdx];

  if (intervalH > oldH || (oldH == 1.0 && intervalH > solver_consts::zeroHLevel))
  {
    mHEstimations[i->problemIdx] = intervalH;
    mNeeRefillQueue = true;
  }
}

template <class FType>
double GOSolver<FType>::CalculateR(const Interval* i)
{
  unsigned problemIdx = i->problemIdx;
  double h = mHEstimations[problemIdx];
  double r = mParameters.r;
  double value = i->delta +
    (i->zr - i->zl) * (i->zr - i->zl) / (i->delta * h * h * r * r) -
      2 * (i->zr + i->zl - 2 * mOptimumEstimations[problemIdx].z) / (r * h);
  return value;
}

template <class FType>
void GOSolver<FType>::InsertIntervals()
{
  for(size_t i = 0; i < mParameters.numThreads; i++)
  {
    //create new interval
    Interval* pNewInterval = new Interval(mNextPoints[i].x, mNextIntervals[i]->xr);
    pNewInterval->zl = mNextPoints[i].z;
    pNewInterval->zr = mNextIntervals[i]->zr;
    pNewInterval->delta = pow(pNewInterval->xr - pNewInterval->xl,
      1. / mProblems.GetDimension());
    pNewInterval->problemIdx = mNextIntervals[i]->problemIdx;
    bool wasInserted =
      mSearchInformations[mNextIntervals[i]->problemIdx].insert(pNewInterval).second;
    if(!wasInserted)
      throw std::runtime_error("Error during interval insertion.");

    //update old interval
    mNextIntervals[i]->xr = mNextPoints[i].x;
    mNextIntervals[i]->zr = mNextPoints[i].z;
    mNextIntervals[i]->delta = pow(mNextIntervals[i]->xr - mNextIntervals[i]->xl,
      1. / mProblems.GetDimension());

    UpdateH(mNextIntervals[i]);
    UpdateH(pNewInterval);

    if(!mNeeRefillQueue)
    {
      pNewInterval->R = CalculateR(pNewInterval);
      mNextIntervals[i]->R = CalculateR(mNextIntervals[i]);
      mQueue.push(pNewInterval);
      mQueue.push(mNextIntervals[i]);
    }
  }
}

template <class FType>
void GOSolver<FType>::MakeTrials()
{
  mNumberOfTrials += mParameters.numThreads;
#pragma omp parallel for num_threads(mParameters.numThreads)
  for(int i = 0; i < (int)mParameters.numThreads; i++)
  {
    mNextPoints[i].z = mProblems.CalculateObjective(mNextPoints[i].y, mNextIntervals[i]->problemIdx);
  }
}

template <class FType>
void GOSolver<FType>::InitDataStructures()
{
  double leftDomainBound[solverMaxDim], rightDomainBound[solverMaxDim];
  mProblems.GetBounds(leftDomainBound, rightDomainBound);
  mEvolvent = Evolvent(mProblems.GetDimension(), mParameters.evloventTightness,
    leftDomainBound, rightDomainBound);
  mQueue = {};
  mNextPoints.resize(mParameters.numThreads);
  mNextIntervals.resize(mParameters.numThreads);
  mActiveProblemsMask.resize(mProblems.GetNumberOfProblems());
  mSearchInformations.resize(mProblems.GetNumberOfProblems());
  std::fill(mActiveProblemsMask.begin(), mActiveProblemsMask.end(), true);
  mHEstimations.resize(mProblems.GetNumberOfProblems());
  std::fill(mHEstimations.begin(), mHEstimations.end(), 1.0);
  mOptimumEstimations.resize(mProblems.GetNumberOfProblems());
  std::fill(mOptimumEstimations.begin(), mOptimumEstimations.end(), Trial(0, HUGE_VAL));
}

template <class FType>
void GOSolver<FType>::ClearDataStructures()
{
  for (size_t i = 0; i < mProblems.GetNumberOfProblems(); i++)
  {
    for(auto it = mSearchInformations[i].begin(); it != mSearchInformations[i].end(); ++it)
      delete *it;
    mSearchInformations[i].clear();
  }
  mQueue = {};
}

template <class FType>
void GOSolver<FType>::FirstIteration()
{
  for (size_t i = 0; i < mProblems.GetNumberOfProblems(); i++)
  {
    Interval* pFirstInterval = new Interval(0., 1.);
    pFirstInterval->delta = 1.;
    pFirstInterval->problemIdx = i;
    double y[solverMaxDim];
    mEvolvent.GetImage(pFirstInterval->xl, y);
    pFirstInterval->zl = mProblems.CalculateObjective(y, i);
    mEvolvent.GetImage(pFirstInterval->xr, y);
    pFirstInterval->zr = mProblems.CalculateObjective(y, i);
    mSearchInformations[i].insert(pFirstInterval);

    if(i < mParameters.numThreads)
    {
      mNextPoints[i].x = 0.5;
      mNextIntervals[i] = pFirstInterval;
    }
  }

  mIterationsCounter = 1;
  mNeeRefillQueue = true;
  mNumberOfActiveProblems = mProblems.GetNumberOfProblems();
  mNumberOfTrials = mProblems.GetNumberOfProblems() * 2;
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
