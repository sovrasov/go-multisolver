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

  void AddProblem(std::shared_ptr<FType> problem)
  {
    mProblems.push_back(problem);
  }

  void DeleteProblem(unsigned problemIndex)
  {
    mProblems.erase(mProblems.begin() + problemIndex);
  }

  void GetBounds(double* lb, double* ub)
  {
    if (mProblems.size() > 0)
      mProblems[0]->GetDomainBounds(lb, ub);
  }

  unsigned Size() const
  {
    return mProblems.size();
  }

  unsigned GetDimension() const
  {
    if(mProblems.size() > 0)
      return mProblems[0]->GetDimension();
    else
      return 0;
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
