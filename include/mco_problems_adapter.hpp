#pragma once

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <exception>


class MultiObjectiveProblemAdapter
{
public:
  using FuncPtr = std::function<double(const double*)>;
protected:
  std::vector<FuncPtr> mObjectives;
  std::vector<FuncPtr> mConstraints;
  std::vector<double> mLBound;
  std::vector<double> mUBound;
  std::vector<double> mLambdas;
  size_t mNumberOfParetoPoints;

public:
  MultiObjectiveProblemAdapter() : mNumberOfParetoPoints(0) {}
  MultiObjectiveProblemAdapter(const std::vector<FuncPtr>& objectives,
                        const std::vector<FuncPtr>& constraints,
                        const std::vector<double>& lb, const std::vector<double>& ub,
                        size_t numParetoPoints = 100) :
                        mObjectives(objectives), mConstraints(constraints),
                        mLBound(lb), mUBound(ub)
  {
    if (mLBound.size() != mUBound.size() || mLBound.size() < 1)
      throw std::runtime_error("Bounds dimension is incorrect");
    if (numParetoPoints < 1)
      throw std::runtime_error("There sould be more than one Pareto point");
    mNumberOfParetoPoints = numParetoPoints;
    if (mNumberOfParetoPoints > 1)
    {
      double h = 1. / (mNumberOfParetoPoints - 1);
      for (size_t i = 0; i < mNumberOfParetoPoints; i++)
        mLambdas.push_back(i*h);
    }
  }

  void SetLambdas(const std::vector<double> lambdas)
  {
    mLambdas = lambdas;
  }

  void GetBounds(double* lb, double* ub, unsigned problemIndex)
  {
    std::copy_n(mLBound.begin(), mLBound.size(), lb);
    std::copy_n(mUBound.begin(), mUBound.size(), ub);
  }

  unsigned GetConstraintsNumber(unsigned problemIndex) const
  {
    return mConstraints.size();
  }

  unsigned GetDimension(unsigned problemIndex = 0) const
  {
    return mUBound.size();
  }

  unsigned GetSize() const
  {
    return mNumberOfParetoPoints;
  }

  double CalculateObjective(const double* y, unsigned problemIdx, unsigned fIndex=0)
  {
    if (fIndex == GetConstraintsNumber(problemIdx))
      return fmax(mObjectives[0](y) * mLambdas[problemIdx], mObjectives[1](y) * (1. - mLambdas[problemIdx]));
    else
      return mConstraints[fIndex](y);
  }

  double GetOptimalValue(unsigned problemIndex) const
  {
    return 0;
  }

  void GetOptimumCoordinates(double* y, unsigned problemIndex) const
  {
  }
};
