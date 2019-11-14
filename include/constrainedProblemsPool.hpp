#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <exception>

template <class FType>
class GCGenProblemsPool
{
protected:

  std::vector<std::shared_ptr<FType>> mProblems;
  std::vector<std::pair<std::vector<double>,std::vector<double>>> mBounds;

public:

  void Add(std::shared_ptr<FType> problem, const std::vector<double>& lBounds,
           const std::vector<double>& uBounds)
  {
    mProblems.push_back(problem);
    if (problem->GetDimension() != lBounds.size() || problem->GetDimension() != uBounds.size())
      throw std::runtime_error("Problem and bound dimension mismatch");
    mBounds.push_back(std::make_pair(lBounds, uBounds));
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
    std::vector<double> tmp_y(y, y + mProblems[problemIndex]->GetDimension());
    if (fIndex == mProblems[problemIndex]->GetConstraintsNumber())
      return mProblems[problemIndex]->ComputeFunction(tmp_y);
    else
      return mProblems[problemIndex]->ComputeConstraint(fIndex, tmp_y);
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
