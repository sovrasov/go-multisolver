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
#include <cmath>

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
  IntervalsQueue mQueue;
  std::vector<std::set<Interval*>> mSearchInformations;

  void InitDataStructures();
  void FirstIteration();
  bool UpdateHConst(Interval*);
  void ClearDataStructures();
	void MakeTrials();
	double CalculateR(Interval*);
	void InsertIntervals();

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

  do {
		MakeTrials();
		InsertIntervals();
		//UpdateHConsts() flag may be setted to true
		//EstimateOptimums() flag may be setted to true
		//Push to queue or refill
		//CalculateNextPoints
		//CheckStop
		mIterationsCounter++;
  } while(false);

  ClearDataStructures();
}

template <class FType>
double GOSolver<FType>::CalculateR(Interval*)
{
	return 0.;
}


template <class FType>
void GOSolver<FType>::InsertIntervals()
{
	for(size_t i = 0; i < mParameters.numThreads; i++)
	{
		//create new interval
		Interval* pNewInterval(mNextPoints[i].x, mNextIntervals[i]->xr);
		pNewInterval->zl = mNextPoints[i].z;
		pNewInterval->zr = mNextIntervals[i]->zr;
		pNewInterval->R = CalculateR(pNewInterval);
		pNewInterval->delta = pow(pNewInterval->xr - pNewInterval->xl,
			1. / mProblems.GetDimension());

		//update old interval
		mNextIntervals[i]->xr = mNextPoints[i].x;
		mNextIntervals[i]->zr = mNextPoints[i].z;
		mNextIntervals[i]->R = CalculateR(mNextIntervals[i]);
		mNextIntervals[i]->delta = pow(mNextIntervals[i]->xr - mNextIntervals[i]->xl,
			1. / mProblems.GetDimension());
		
		//updatehCponst
	}
}

template <class FType>
void GOSolver<FType>::MakeTrials()
{
//#pragma omp parallel for num_threads(mParameters.numThreads)
	for(size_t i = 0; i < mParameters.numThreads; i++)
	{
		mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y);
		mNextPoints[i].z = mProblems.CalculateObjective(mNextPoints[i].y, mNextIntervals[i]->problemIdx);
	}
}

template <class FType>
void GOSolver<FType>::InitDataStructures()
{
  mQueue.Clear();
  mEvolvent = Evolvent(mProblems.GetDimension(), mParameters.evloventTightness);
  mOptimumEstimations.resize(mProblems.GetNumberOfProblems());
  mNextPoints.resize(mParameters.numThreads);
  mNextIntervals.resize(mParameters.numThreads);
  mActiveProblemsMask.resize(mProblems.GetNumberOfProblems());
  mSearchInformations.resize(mProblems.GetNumberOfProblems());
  std::fill(mActiveProblemsMask.begin(), mActiveProblemsMask.end(), true);
	mHEstimations.resize(mProblems.GetNumberOfProblems());
	std::fill(mHEstimations.begin(), mHEstimations.end(), 1.0);
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
  mQueue.Clear();
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
