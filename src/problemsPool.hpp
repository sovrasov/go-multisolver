#ifndef PROBLEMS_POOL_HPP
#define PROBLEMS_POOL_HPP

#include <vector>

template <class FType>
class ProblemsPool
{
protected:
  std::vector<FType*> mProblems;

public:

  void AddProblem(FType* problem)
  {
    mProblems.push_back(problem);
  }

  void GetBounds(double* lb, double* ub)
  {
    if (mProblems.size() > 0)
      mProblems[0]->GetBounds(lb, ub);
  }

  unsigned GetNumberOfProblems() const
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
};

#endif
