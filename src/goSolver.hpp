#ifndef GO_SOLVER_HPP
#define GO_SOLVER_HPP

#include "problemsPool.hpp"
#include "dataTypes.hpp"
#include "evolvent.hpp"

#include <vector>
#include <set>
#include <queue>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace solver_internal
{
  const double zeroHLevel = 1e-12;
  using PriorityQueue = std::priority_queue<Interval*, std::vector<Interval*>, CompareByR>;

  bool checkVectorsDiff(const double* y1, const double* y2, size_t dim, double eps)
  {
    for (size_t i = 0; i < dim; i++)
    {
      if (fabs(y1[i] - y2[i]) > eps)
        return true;
    }

    return false;
  }

  double vectorsMaxDiff(const double* y1, const double* y2, size_t dim)
  {
    double diff = 0.;
    for (size_t i = 0; i < dim; i++)
      diff = fmax(fabs(y1[i] - y2[i]), diff);

    return diff;
  }
}

enum class StopType {Accuracy, OptimumVicinity, OptimalValue};

struct SolverParameters
{
  double eps;
  double r;
  unsigned numThreads;
  unsigned iterationsLimit;
  unsigned evloventTightness = 12;
  StopType criterion;
  bool logDeviations;

  SolverParameters() {}
  SolverParameters(double _eps, double _r,
      unsigned _numThreads, unsigned _iterationsLimit, StopType _criterion = StopType::Accuracy) :
        eps(_eps), r(_r), numThreads(_numThreads), iterationsLimit(_iterationsLimit),
        criterion(_criterion)
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
  std::vector<Trial> mOptimumEstimations;
  std::vector<Interval*> mNextIntervals;
  std::vector<Trial> mNextPoints;
  std::vector<std::set<Interval*>> mSearchInformations;
  solver_internal::PriorityQueue mQueue;
  std::vector<StatPoint> mStatiscics;
  Evolvent mEvolvent;
  bool mNeeRefillQueue;
  unsigned mIterationsCounter;
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
  void CollectStatistics();

public:

  void SetProblemsPool(ProblemsPool<FType>& problems);
  void SetParameters(SolverParameters& params);
  void Solve();
  std::vector<Trial> GetOptimumEstimations();
  unsigned GetTrialsNumber() const { return mNumberOfTrials; }
  unsigned GetIterationsNumber() const { return mIterationsCounter; }
  std::vector<StatPoint> GetStatistics() const { return mStatiscics; }
};

template <class FType>
void GOSolver<FType>::Solve()
{
  bool needStop = false;
  InitDataStructures();
  FirstIteration();

  do {
    MakeTrials();
    EstimateOptimums();
    InsertIntervals();
    if(mParameters.logDeviations && mNumberOfTrials % 1000 == 0)
      CollectStatistics();
    if (mNeeRefillQueue || mQueue.size() < mParameters.numThreads)
      RefillQueue();
    CalculateNextPoints();
    needStop = CheckStopCondition();
    mIterationsCounter++;
  } while(mIterationsCounter < mParameters.iterationsLimit && !needStop);

  MakeTrials();
  EstimateOptimums();
  if(mParameters.logDeviations)
    CollectStatistics();

  ClearDataStructures();
}

template <class FType>
void GOSolver<FType>::CollectStatistics()
{
  StatPoint currentDevs(mNumberOfTrials, 0., 0.);
  for (size_t j = 0; j < mProblems.Size(); j++)
  {
    double optPoint[solverMaxDim];
    mProblems.GetOptimumCoordinates(optPoint, j);
    double difference = solver_internal::vectorsMaxDiff(optPoint, mOptimumEstimations[j].y, mProblems.GetDimension());
    //if (mSearchInformations[j].size() == 1)
      //difference = ;
    currentDevs.meanDev += difference;
    currentDevs.maxDev = fmax(currentDevs.maxDev, difference);
  }
  currentDevs.meanDev /= mProblems.Size();

  mStatiscics.push_back(currentDevs);
}


template <class FType>
bool GOSolver<FType>::CheckStopCondition()
{
  bool needStop = false;
  for (size_t i = 0; i < mParameters.numThreads; i++)
  {
    bool isOptimumReached = false;
    switch (mParameters.criterion)
    {
    case StopType::Accuracy:
      isOptimumReached = mNextIntervals[i]->delta < mParameters.eps;
      break;
    case StopType::OptimumVicinity:
    {
      double optimum[solverMaxDim];
      mProblems.GetOptimumCoordinates(optimum, mNextIntervals[i]->problemIdx);
      isOptimumReached = !solver_internal::checkVectorsDiff(
        optimum, mNextPoints[i].y, mProblems.GetDimension(), mParameters.eps);
    }
      break;
    case StopType::OptimalValue:
      isOptimumReached = mOptimumEstimations[mNextIntervals[i]->problemIdx].z -
        mProblems.GetOptimalValue(mNextIntervals[i]->problemIdx) < mParameters.eps;
    }

    if (isOptimumReached)
    {
      mActiveProblemsMask[mNextIntervals[i]->problemIdx] = false;
      mNumberOfActiveProblems--;
      std::cout << "Problem # " << mNextIntervals[i]->problemIdx + 1 <<
        " has been solved! Trials performed: " << mNumberOfTrials <<
        " Problems left: " << mNumberOfActiveProblems << "\n";
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
        (((mNextIntervals[i]->zr - mNextIntervals[i]->zl) > 0.) ? 1. : -1.) *
        pow(fabs(mNextIntervals[i]->zr - mNextIntervals[i]->zl) /
        mHEstimations[mNextIntervals[i]->problemIdx], mProblems.GetDimension()) / 2. / mParameters.r;

    if (mNextPoints[i].x > mNextIntervals[i]->xr || mNextPoints[i].x < mNextIntervals[i]->xl)
      throw std::runtime_error("The next point is outside of the subdivided interval");

    mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y);
  }
}

template <class FType>
void GOSolver<FType>::RefillQueue()
{
  mQueue = solver_internal::PriorityQueue();

  for(size_t i = 0; i < mProblems.Size(); i++)
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
  double intervalH = fabs(i->zr - i->zl) / i->delta;
  double oldH = mHEstimations[i->problemIdx];

  if (intervalH > oldH || (oldH == 1.0 && intervalH > solver_internal::zeroHLevel))
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
  mQueue = solver_internal::PriorityQueue();
  mNextPoints.resize(mParameters.numThreads);
  mNextIntervals.resize(mParameters.numThreads);
  mActiveProblemsMask.resize(mProblems.Size());
  mSearchInformations.resize(mProblems.Size());
  mStatiscics.resize(0);
  std::fill(mActiveProblemsMask.begin(), mActiveProblemsMask.end(), true);
  mHEstimations.resize(mProblems.Size());
  std::fill(mHEstimations.begin(), mHEstimations.end(), 1.0);
  mOptimumEstimations.resize(mProblems.Size());
  std::fill(mOptimumEstimations.begin(), mOptimumEstimations.end(), Trial(0., HUGE_VAL));
}

template <class FType>
void GOSolver<FType>::ClearDataStructures()
{
  for (size_t i = 0; i < mProblems.Size(); i++)
  {
    for(auto it = mSearchInformations[i].begin(); it != mSearchInformations[i].end(); ++it)
      delete *it;
    mSearchInformations[i].clear();
  }
  mQueue = solver_internal::PriorityQueue();
}

template <class FType>
void GOSolver<FType>::FirstIteration()
{
  for (size_t i = 0; i < mProblems.Size(); i++)
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
    UpdateH(pFirstInterval);

    if(i < mParameters.numThreads)
    {
      mNextPoints[i].x = 0.5;
      mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y);
      mNextIntervals[i] = pFirstInterval;
    }
  }

  mIterationsCounter = 1;
  mNeeRefillQueue = true;
  mNumberOfActiveProblems = mProblems.Size();
  mNumberOfTrials = mProblems.Size() * 2;
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
