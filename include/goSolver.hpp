#ifndef GO_SOLVER_HPP
#define GO_SOLVER_HPP

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
  double epsR = 0;
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

template <class PoolType>
class GOSolver
{
protected:

  SolverParameters mParameters;
  PoolType mProblems;
  std::vector<bool> mActiveProblemsMask;
  std::vector<std::vector<double>> mHEstimations;
  std::vector<std::vector<double>> mZEstimations;
  std::vector<Trial> mOptimumEstimations;
  std::vector<Interval*> mNextIntervals;
  std::vector<Trial> mNextPoints;
  std::vector<double> mMinDifferences;
  std::vector<int> mMaxIndexes;
  std::vector<std::set<Interval*>> mSearchInformations;
  std::vector<std::vector<double>> mLowerDomainBounds;
  std::vector<std::vector<double>> mUpperDomainBounds;
  solver_internal::PriorityQueue mQueue;
  std::vector<StatPoint> mStatiscics;
  Evolvent mEvolvent;
  bool mNeedRefillQueue;
  unsigned mIterationsCounter;
  unsigned mNumberOfActiveProblems;
  unsigned mNumberOfTrials;
  double mDimExponent;

  void InitDataStructures();
  void FirstIteration();
  void ClearDataStructures();
  void MakeTrials();
  void MakeTrial(Trial&, int);
  void MakeTrial(Point&, const double*, int);
  void MakeTrial(const double*, int, double*, int&);
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

  void SetProblemsPool(PoolType& problems);
  void SetParameters(SolverParameters& params);
  void Solve();
  std::vector<Trial> GetOptimumEstimations();
  unsigned GetTrialsNumber() const { return mNumberOfTrials; }
  unsigned GetIterationsNumber() const { return mIterationsCounter; }
  std::vector<StatPoint> GetStatistics() const { return mStatiscics; }
};

template <class PoolType>
double GOSolver<PoolType>::GetNextPointCoordinate(const Interval* i) const
{
  int v = i->xl.v;
  double diff = i->xr.z[v] - i->xl.z[v];
  return 0.5 * (i->xl.x + i->xr.x) -
    ((diff > 0.) ? 1. : -1.) * pow(fabs(diff) /
      mHEstimations[i->problemIdx][v], mProblems.GetDimension()) / 2. / mParameters.r;
}

template <class PoolType>
void GOSolver<PoolType>::Solve()
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
    if (mNeedRefillQueue || mQueue.size() < mParameters.numThreads)
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

template <class PoolType>
void GOSolver<PoolType>::CollectStatistics()
{
  StatPoint currentDevs(mNumberOfTrials);
  for (size_t j = 0; j < mProblems.GetSize(); j++)
  {
    double difference;
    if (mParameters.criterion != StopType::OptimalValue)
      difference = mMinDifferences[j];
    else
      difference = mOptimumEstimations[j].GetZ() - mProblems.GetOptimalValue(j);
    currentDevs.meanDev += difference;
    currentDevs.maxDev = fmax(currentDevs.maxDev, difference);
    if (difference < mParameters.eps)
      currentDevs.problems_solved++;
  }
  currentDevs.meanDev /= mProblems.GetSize();

  mStatiscics.push_back(currentDevs);
}

template <class PoolType>
bool GOSolver<PoolType>::CheckStopCondition()
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
      isOptimumReached = mNextPoints[i].GetZ() -
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
          mNextPoints[i].GetZ() - mProblems.GetOptimalValue(mNextIntervals[i]->problemIdx) <<
          " H estimation: " << mHEstimations[mNextIntervals[i]->problemIdx][0] << "\n";
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

template <class PoolType>
void GOSolver<PoolType>::CalculateNextPoints()
{
  for(size_t i = 0; i < mParameters.numThreads; i++)
  {
    mNextIntervals[i] = mQueue.top();
    mQueue.pop();
    mNextPoints[i].x = GetNextPointCoordinate(mNextIntervals[i]);

    if (mNextPoints[i].x >= mNextIntervals[i]->xr.x || mNextPoints[i].x <= mNextIntervals[i]->xl.x)
      throw std::runtime_error("The next point is outside of the subdivided interval");

    mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y,
      mLowerDomainBounds[mNextIntervals[i]->problemIdx].data(),
      mUpperDomainBounds[mNextIntervals[i]->problemIdx].data());
  }
}

template <class PoolType>
void GOSolver<PoolType>::RefillQueue()
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

  mNeedRefillQueue = false;
}

template <class PoolType>
void GOSolver<PoolType>::EstimateOptimums()
{
  for(size_t i = 0; i < mParameters.numThreads; i++)
  {
    unsigned problemIdx = mNextIntervals[i]->problemIdx;
    if(mNextPoints[i].z < mOptimumEstimations[problemIdx].z)
    {
      mOptimumEstimations[problemIdx] = mNextPoints[i];
      mNeedRefillQueue = true;
    }
  }
}

template <class PoolType>
void GOSolver<PoolType>::UpdateH(const Interval* i)
{
  double intervalH = fabs(i->xr.z[i->xr.v] - i->xl.z[i->xl.v]) / i->delta;
  double oldH = mHEstimations[i->problemIdx][0];

  if (intervalH > oldH || (oldH == 1.0 && intervalH > solver_internal::zeroHLevel))
  {
    mHEstimations[i->problemIdx][0] = intervalH;
    mNeedRefillQueue = true;
  }
}

template <class PoolType>
double GOSolver<PoolType>::CalculateR(const Interval* i)
{
  unsigned problemIdx = i->problemIdx;
  double h = mHEstimations[problemIdx][0];
  double r = mParameters.r;
  int v = i->xr.v;
  double value = i->delta +
    (i->xr.z[v] - i->xl.z[v]) * (i->xr.z[v] - i->xl.z[v]) / (i->delta * h * h * r * r) -
      2 * (i->xr.z[v] + i ->xl.z[v] - 2 * mOptimumEstimations[problemIdx].z[v]) / (r * h);
  return value;
}

template <class PoolType>
void GOSolver<PoolType>::InsertIntervals()
{
  for(size_t i = 0; i < mParameters.numThreads; i++)
  {
    Interval* pNewInterval = new Interval(mNextPoints[i], mNextIntervals[i]->xr);
    pNewInterval->delta = pow(pNewInterval->xr.x - pNewInterval->xl.x, mDimExponent);
    pNewInterval->problemIdx = mNextIntervals[i]->problemIdx;
    bool wasInserted =
      mSearchInformations[mNextIntervals[i]->problemIdx].insert(pNewInterval).second;
    if(!wasInserted)
      throw std::runtime_error("Error during interval insertion.");

    mNextIntervals[i]->xr = mNextPoints[i];
    mNextIntervals[i]->delta = pow(mNextIntervals[i]->xr.x - mNextIntervals[i]->xl.x,
      mDimExponent);

    UpdateH(mNextIntervals[i]);
    UpdateH(pNewInterval);

    if(!mNeedRefillQueue)
    {
      pNewInterval->R = CalculateR(pNewInterval);
      mNextIntervals[i]->R = CalculateR(mNextIntervals[i]);
      mQueue.push(pNewInterval);
      mQueue.push(mNextIntervals[i]);
    }
  }
}

template <class PoolType>
void GOSolver<PoolType>::MakeTrials()
{
  mNumberOfTrials += mParameters.numThreads;
#pragma omp parallel for num_threads(mParameters.numThreads)
  for(int i = 0; i < (int)mParameters.numThreads; i++)
  {
    MakeTrial(mNextPoints[i], mNextIntervals[i]->problemIdx);
  }
}

template <class PoolType>
void GOSolver<PoolType>::MakeTrial(Trial& t, int problem_idx)
{
  int index;
  MakeTrial(t.y, problem_idx, t.z, index);
  t.v = index;
}

template <class PoolType>
void GOSolver<PoolType>::MakeTrial(Point& p, const double* y, int problem_idx)
{
  int index;
  MakeTrial(y, problem_idx, p.z, index);
  p.v = index;
}

template <class PoolType>
void GOSolver<PoolType>::MakeTrial(const double* y, int problem_idx, double* z, int& idx)
{
  idx = 0;
  while(idx < mProblems.GetConstraintsNumber(problem_idx))
  {
    double val = mProblems.CalculateObjective(y, problem_idx, idx);
    z[idx] = val;
    if (val > 0)
      break;
    idx++;
  }

  if(idx > mMaxIndexes[problem_idx])
  {
    mMaxIndexes[problem_idx] = idx;
    for(int i = 0; i < mMaxIndexes[problem_idx]; i++)
      mZEstimations[problem_idx][i] = -mParameters.epsR*mHEstimations[problem_idx][i];
    mNeedRefillQueue = true;
  }

  if(idx == mProblems.GetConstraintsNumber(problem_idx))
    z[idx] = mProblems.CalculateObjective(y, problem_idx, idx);

  if(idx == mMaxIndexes[problem_idx] &&
     z[idx] < mZEstimations[problem_idx][idx])
  {
    mZEstimations[problem_idx][idx] = z[idx];
    mNeedRefillQueue = true;
  }
}

template <class PoolType>
void GOSolver<PoolType>::InitDataStructures()
{
  double leftDomainBound[solverMaxDim], rightDomainBound[solverMaxDim];
  mProblems.GetBounds(leftDomainBound, rightDomainBound, 0);
  mEvolvent = Evolvent(mProblems.GetDimension(), mParameters.evloventTightness);

  mQueue = solver_internal::PriorityQueue();
  mNextPoints.resize(mParameters.numThreads);
  mNextIntervals.resize(mParameters.numThreads);
  mActiveProblemsMask.resize(mProblems.GetSize());
  mSearchInformations.resize(mProblems.GetSize());
  mLowerDomainBounds.resize(mProblems.GetSize());
  mUpperDomainBounds.resize(mProblems.GetSize());
  mMaxIndexes.resize(mProblems.GetSize());
  std::fill(mMaxIndexes.begin(), mMaxIndexes.end(), -1);
  mStatiscics.resize(0);
  std::fill(mActiveProblemsMask.begin(), mActiveProblemsMask.end(), true);
  mHEstimations.resize(mProblems.GetSize());
  std::fill(mHEstimations.begin(), mHEstimations.end(),
            std::vector<double>(solverMaxFunctionsNum, 1.0));
  mZEstimations.resize(mProblems.GetSize());
  std::fill(mZEstimations.begin(), mZEstimations.end(),
            std::vector<double>(solverMaxFunctionsNum, std::numeric_limits<double>::max()));
  mOptimumEstimations.resize(mProblems.GetSize());
  std::fill(mOptimumEstimations.begin(), mOptimumEstimations.end(), Trial(0.));
  mMinDifferences.resize(mProblems.GetSize());
  std::fill(mMinDifferences.begin(), mMinDifferences.end(), solver_internal::
    vectorsMaxDiff(leftDomainBound, rightDomainBound, mProblems.GetDimension()));

  mDimExponent = 1. / mProblems.GetDimension();
}

template <class PoolType>
void GOSolver<PoolType>::ClearDataStructures()
{
  for (size_t i = 0; i < mProblems.GetSize(); i++)
  {
    for(auto it = mSearchInformations[i].begin(); it != mSearchInformations[i].end(); ++it)
      delete *it;
    mSearchInformations[i].clear();
  }
  mQueue = solver_internal::PriorityQueue();
}

template <class PoolType>
void GOSolver<PoolType>::FirstIteration()
{
  for (size_t i = 0; i < mProblems.GetSize(); i++)
  {
    mLowerDomainBounds[i].resize(mProblems.GetDimension());
    mUpperDomainBounds[i].resize(mProblems.GetDimension());
    mProblems.GetBounds(mLowerDomainBounds[i].data(), mUpperDomainBounds[i].data(), i);

    Interval* pFirstInterval = new Interval(0., 1.);
    pFirstInterval->delta = 1.;
    pFirstInterval->problemIdx = i;
    double yl[solverMaxDim], yr[solverMaxDim];
    mEvolvent.GetImage(pFirstInterval->xl.x, yl, mLowerDomainBounds[i].data(), mUpperDomainBounds[i].data());
    MakeTrial(pFirstInterval->xl, yl, i);
    mEvolvent.GetImage(pFirstInterval->xr.x, yr, mLowerDomainBounds[i].data(), mUpperDomainBounds[i].data());
    MakeTrial(pFirstInterval->xr, yr, i);
    mSearchInformations[i].insert(pFirstInterval);
    UpdateH(pFirstInterval);
    if(pFirstInterval->xl.z[pFirstInterval->xl.v] < pFirstInterval->xr.z[pFirstInterval->xr.v] &&
       pFirstInterval->xl.v >= pFirstInterval->xr.v)
    {
        mOptimumEstimations[i] = Trial(pFirstInterval->xl.x);
        std::copy_n(pFirstInterval->xl.z, mProblems.GetDimension(), mOptimumEstimations[i].z);
        std::copy_n(yl, mProblems.GetDimension(), mOptimumEstimations[i].y);
    }
    else
    {
      mOptimumEstimations[i] = Trial(pFirstInterval->xr.x);
      std::copy_n(pFirstInterval->xr.z, mProblems.GetDimension(), mOptimumEstimations[i].z);
      std::copy_n(yr, mProblems.GetDimension(), mOptimumEstimations[i].y);
    }
  }

  RefillQueue();
  CalculateNextPoints();

  mIterationsCounter = 1;
  mNeedRefillQueue = true;
  mNumberOfActiveProblems = mProblems.GetSize();
  mNumberOfTrials = mProblems.GetSize() * 2;
}

template <class PoolType>
void GOSolver<PoolType>::SetProblemsPool(PoolType& problems)
{
  mProblems = problems;
}

template <class PoolType>
void GOSolver<PoolType>::SetParameters(SolverParameters& params)
{
  mParameters = params;
}

template <class PoolType>
std::vector<Trial> GOSolver<PoolType>::GetOptimumEstimations()
{
  return mOptimumEstimations;
}

#endif
