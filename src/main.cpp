#include "Grishagin/grishagin_function.hpp"
#include "GKLS/gkls_function.hpp"
#include "problemsPool.hpp"
#include "goSolver.hpp"

#include <iostream>

int main(int argc, const char** argv)
{
  ProblemsPool<gkls::GKLSFunction> pool;
  gkls::GKLSFunction** functions = new gkls::GKLSFunction*[100];
  for(unsigned i = 0; i < 1; i++)
  {
    gkls::GKLSFunction* func = new gkls::GKLSFunction();
    func->SetFunctionClass(gkls::Simple, 2);
    func->SetType(gkls::TD);
    func->SetFunctionNumber(i + 1);
    pool.AddProblem(func);
  }

  std::cout << "Problems pool created\n";

  SolverParameters parameters(0.01, 4.0, 1, 100000);
  GOSolver<gkls::GKLSFunction> solver;
  solver.SetParameters(parameters);
  solver.SetProblemsPool(pool);
  std::cout << "Solver started\n";
  solver.Solve();
  std::vector<Trial> optimumEstimations = solver.GetOptimumEstimations();

	for(size_t i = 0; i < optimumEstimations.size(); i++)
		std::cout << "Optimum value in problem #" << i + 1 << ": " << optimumEstimations[i].z << "\n";

  for(int i = 0; i < 100; i++)
    delete functions[i];
  delete[] functions;

  return 0;
}
