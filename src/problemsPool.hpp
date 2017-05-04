#ifndef PROBLEMS_POOL_HPP
#define PROBLEMS_POOL_HPP

#include <vector>
#include <memory>

template <class FType>
class ProblemsPool
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
    mProblems[problemIndex]->GetDomainBounds(lb, ub);
  }

  unsigned GetSize() const
  {
    return mProblems.size();
  }

  unsigned GetDimension(unsigned problemIndex = 0) const
  {
    return mProblems[problemIndex]->GetDimension();
  }

  double CalculateObjective(const double* y, unsigned problemIndex)
  {
    return mProblems[problemIndex]->Calculate(y);
  }

  double GetOptimalValue(unsigned problemIndex) const
  {
    return mProblems[problemIndex]->GetOptimalValue();
  }

  void GetOptimumCoordinates(double* y, unsigned problemIndex) const
  {
    mProblems[problemIndex]->GetOptimumCoordinates(y);
  }
};

#endif
