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
  unsigned trialsLimit;
  unsigned evloventTightness = 12;
  StopType criterion;
  bool logDeviations;
  bool verbose = false;
  unsigned statisticsUpdateStep = 1000;

  SolverParameters() {}
  SolverParameters(double _eps, double _r,
      unsigned _numThreads, unsigned _trialsLimit, StopType _criterion = StopType::Accuracy) :
        eps(_eps), r(_r), numThreads(_numThreads), trialsLimit(_trialsLimit),
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
  std::vector<double> mMinDifferences;
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
  double GetNextPointCoordinate(const Interval*) const;

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
double GOSolver<FType>::GetNextPointCoordinate(const Interval* i) const
{
  return 0.5 * (i->xl + i->xr) -
    (((i->zr - i->zl) > 0.) ? 1. : -1.) * pow(fabs(i->zr - i->zl) /
      mHEstimations[i->problemIdx], mProblems.GetDimension()) / 2. / mParameters.r;
}

template <class FType>
void GOSolver<FType>::Solve()
{
  bool needStop = false;
  InitDataStructures();
  FirstIteration();
  MakeTrials();
  if(mParameters.logDeviations)
    CollectStatistics();

  do {
    EstimateOptimums();
    InsertIntervals();
    if(mParameters.logDeviations &&
        mNumberOfTrials % mParameters.statisticsUpdateStep == 0)
      CollectStatistics();
    if (mNeeRefillQueue || mQueue.size() < mParameters.numThreads)
      RefillQueue();
    CalculateNextPoints();
    MakeTrials();
    needStop = CheckStopCondition();
    mIterationsCounter++;
  } while(mNumberOfTrials < mParameters.trialsLimit && !needStop);

  if(mParameters.logDeviations)
    CollectStatistics();

  ClearDataStructures();
}

template <class FType>
void GOSolver<FType>::CollectStatistics()
{
  StatPoint currentDevs(mNumberOfTrials, 0., 0.);
  for (size_t j = 0; j < mProblems.GetSize(); j++)
  {
    double difference;
    if (mParameters.criterion != StopType::OptimalValue)
      difference = mMinDifferences[j];
    else
      difference = mOptimumEstimations[j].z - mProblems.GetOptimalValue(j);
    currentDevs.meanDev += difference;
    currentDevs.maxDev = fmax(currentDevs.maxDev, difference);
  }
  currentDevs.meanDev /= mProblems.GetSize();

  mStatiscics.push_back(currentDevs);
}

template <class FType>
bool GOSolver<FType>::CheckStopCondition()
{
  bool needStop = false;
  bool needRenewIntervals = false;
  for (size_t i = 0; i < mParameters.numThreads && mNumberOfActiveProblems > 0; i++)
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
      double optDistance = solver_internal::vectorsMaxDiff(optimum, mNextPoints[i].y,
        mProblems.GetDimension());
      isOptimumReached = optDistance < mParameters.eps;
      if(mParameters.logDeviations)
        mMinDifferences[mNextIntervals[i]->problemIdx] = std::min(optDistance,
          mMinDifferences[mNextIntervals[i]->problemIdx]);
    }
      break;
    case StopType::OptimalValue:
      isOptimumReached = mNextPoints[i].z -
        mProblems.GetOptimalValue(mNextIntervals[i]->problemIdx) < mParameters.eps;
    }

    if (isOptimumReached && mActiveProblemsMask[mNextIntervals[i]->problemIdx])
    {
      mActiveProblemsMask[mNextIntervals[i]->problemIdx] = false;
      mNumberOfActiveProblems--;

      mOptimumEstimations[mNextIntervals[i]->problemIdx] = mNextPoints[i];
      //TODO: use all search data to estimate optimum

      if (mParameters.verbose)
      {
        double optimum[solverMaxDim];
        mProblems.GetOptimumCoordinates(optimum, mNextIntervals[i]->problemIdx);
        std::cout << "Problem # " << mNextIntervals[i]->problemIdx + 1 <<
          " has been solved! Trials performed: " << mNumberOfTrials <<
          " Problems left: " << mNumberOfActiveProblems << " coordinates diff: " <<
          solver_internal::vectorsMaxDiff(
            optimum, mNextPoints[i].y, mProblems.GetDimension()) << " Values diff: " <<
          mNextPoints[i].z - mProblems.GetOptimalValue(mNextIntervals[i]->problemIdx) <<
          " H estimation: " << mHEstimations[mNextIntervals[i]->problemIdx] << "\n";
      }

      if (mNumberOfActiveProblems)
        needRenewIntervals = true;
    }
  }

  if (needRenewIntervals)
  {
    RefillQueue();
    CalculateNextPoints();
    MakeTrials();
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
    mNextPoints[i].x = GetNextPointCoordinate(mNextIntervals[i]);

    if (mNextPoints[i].x >= mNextIntervals[i]->xr || mNextPoints[i].x <= mNextIntervals[i]->xl)
      throw std::runtime_error("The next point is outside of the subdivided interval");

    mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y);
  }
}

template <class FType>
void GOSolver<FType>::RefillQueue()
{
  mQueue = solver_internal::PriorityQueue();

  for(size_t i = 0; i < mProblems.GetSize(); i++)
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
  mActiveProblemsMask.resize(mProblems.GetSize());
  mSearchInformations.resize(mProblems.GetSize());
  mStatiscics.resize(0);
  std::fill(mActiveProblemsMask.begin(), mActiveProblemsMask.end(), true);
  mHEstimations.resize(mProblems.GetSize());
  std::fill(mHEstimations.begin(), mHEstimations.end(), 1.0);
  mOptimumEstimations.resize(mProblems.GetSize());
  std::fill(mOptimumEstimations.begin(), mOptimumEstimations.end(), Trial(0., HUGE_VAL));
  mMinDifferences.resize(mProblems.GetSize());
  std::fill(mMinDifferences.begin(), mMinDifferences.end(), solver_internal::
    vectorsMaxDiff(leftDomainBound, rightDomainBound, mProblems.GetDimension()));
}

template <class FType>
void GOSolver<FType>::ClearDataStructures()
{
  for (size_t i = 0; i < mProblems.GetSize(); i++)
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
  for (size_t i = 0; i < mProblems.GetSize(); i++)
  {
    Interval* pFirstInterval = new Interval(0., 1.);
    pFirstInterval->delta = 1.;
    pFirstInterval->problemIdx = i;
    double yl[solverMaxDim], yr[solverMaxDim];
    mEvolvent.GetImage(pFirstInterval->xl, yl);
    pFirstInterval->zl = mProblems.CalculateObjective(yl, i);
    mEvolvent.GetImage(pFirstInterval->xr, yr);
    pFirstInterval->zr = mProblems.CalculateObjective(yr, i);
    mSearchInformations[i].insert(pFirstInterval);
    UpdateH(pFirstInterval);
    if(pFirstInterval->zl < pFirstInterval->zr)
    {
      mOptimumEstimations[i] = Trial(pFirstInterval->xl, pFirstInterval->zl);
      std::copy_n(yl, mProblems.GetDimension(), mOptimumEstimations[i].y);
    }
    else
    {
      mOptimumEstimations[i] = Trial(pFirstInterval->xr, pFirstInterval->zr);
      std::copy_n(yr, mProblems.GetDimension(), mOptimumEstimations[i].y);
    }
  }

  RefillQueue();
  CalculateNextPoints();

  mIterationsCounter = 1;
  mNeeRefillQueue = true;
  mNumberOfActiveProblems = mProblems.GetSize();
  mNumberOfTrials = mProblems.GetSize() * 2;
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
