#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <exception>
#include <functional>

template <class FType>
class GCGenProblemsPool
{
protected:

  std::vector<std::shared_ptr<FType>> mProblems;
  std::vector<std::pair<std::vector<double>,std::vector<double>>> mBounds;
  std::function<double()> mComputeLoad = [] { return 1; };

public:

  void Add(std::shared_ptr<FType> problem)
  {
    mProblems.push_back(problem);
    std::vector<double> lb, ub;
    problem->GetBounds(lb, ub);
    if (problem->GetDimension() != lb.size() || problem->GetDimension() != ub.size())
      throw std::runtime_error("Problem and bound dimension mismatch");
    mBounds.push_back(std::make_pair(lb, ub));
  }

  void SetComputeLoad(std::function<double()> compute)
  {
    mComputeLoad = compute;
  }

  void GetBounds(double* lb, double* ub, unsigned problemIndex)
  {
    std::copy_n(mBounds[problemIndex].first.begin(), mBounds[problemIndex].first.size(), lb);
    std::copy_n(mBounds[problemIndex].second.begin(), mBounds[problemIndex].second.size(), ub);
  }

  unsigned GetConstraintsNumber(unsigned problemIndex) const
  {
    return static_cast<size_t>(mProblems[problemIndex]->GetConstraintsNumber());
  }

  unsigned GetDimension(unsigned problemIndex = 0) const
  {
    return mProblems[problemIndex]->GetDimension();
  }

  unsigned GetSize() const
  {
    return mProblems.size();
  }

  double CalculateObjective(const double* y, unsigned problemIndex, unsigned fIndex=0)
  {
    double k = mComputeLoad();
    std::vector<double> tmp_y(y, y + mProblems[problemIndex]->GetDimension());
    if (fIndex == mProblems[problemIndex]->GetConstraintsNumber())
    {
      double val = k * mProblems[problemIndex]->ComputeFunction(tmp_y);
      return val / k;
    }
    else
    {
      double val = k * mProblems[problemIndex]->ComputeConstraint(fIndex, tmp_y);
      return val / k;
    }
  }

  double GetOptimalValue(unsigned problemIndex) const
  {
    return mProblems[problemIndex]->GetOptimumValue();
  }

  void GetOptimumCoordinates(double* y, unsigned problemIndex) const
  {
    auto point = mProblems[problemIndex]->GetOptimumPoint();
    std::copy(point.data(), point.data() + point.size(), y);
  }
};
