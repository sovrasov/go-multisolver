#pragma once

#include <vector>
#include <memory>

template <class FType>
class GCGenProblemsPool
{
protected:

  std::vector<std::shared_ptr<FType>> mProblems;

public:

  void Add(std::shared_ptr<FType> problem)
  {
    mProblems.push_back(problem);
  }

  void GetBounds(double* lb, double* ub, unsigned problemIndex)
  {
    throw -1;
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
