#include "Grishagin/grishagin_function.hpp"
#include "GKLS/gkls_function.hpp"
#include "problemsPool.hpp"
#include "goSolver.hpp"

#include <iostream>
#include <chrono>
#include <memory>

int main(int argc, const char** argv)
{
  const unsigned nProblems = 100;
  ProblemsPool<gkls::GKLSFunction> pool;
  for(unsigned i = 0; i < nProblems; i++)
  {
    gkls::GKLSFunction* func = new gkls::GKLSFunction();
    func->SetFunctionClass(gkls::Simple, 4);
    func->SetType(gkls::TD);
    func->SetFunctionNumber(i + 1);
    pool.AddProblem(std::shared_ptr<gkls::GKLSFunction>(func));
  }

  std::cout << "Problems pool created\n";

  SolverParameters parameters(0.02, 4.7, 4, 5000000, StopType::OptimumVicinity);
  GOSolver<gkls::GKLSFunction> solver;
  solver.SetParameters(parameters);
  solver.SetProblemsPool(pool);
  std::cout << "Solver started\n";
  auto start = std::chrono::system_clock::now();
  solver.Solve();
  auto end = std::chrono::system_clock::now();
  std::vector<Trial> optimumEstimations = solver.GetOptimumEstimations();

	//for(size_t i = 0; i < optimumEstimations.size(); i++)
	//	std::cout << "Optimum value in problem #" << i + 1 << ": " << optimumEstimations[i].z << "\n";

  std::cout << "Number of trials: " << solver.GetTrialsNumber()
    << "\nNumber of iterations: " << solver.GetIterationsNumber() << "\n";

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Time elapsed: " << elapsed_seconds.count() << "s\n";

  return 0;
}
