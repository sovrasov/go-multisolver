#include "Grishagin/grishagin_function.hpp"
#include "GKLS/gkls_function.hpp"
#include "problemsPool.hpp"
#include "goSolver.hpp"

#include <iostream>
#include <chrono>
#include <memory>
#include <cmdline.h>

int main(int argc, char** argv)
{
  cmdline::parser parser;
  parser.add<int>("dimension", 'd', "test problems dumension", false, 2,
     cmdline::range(1, 5));
  parser.add<int>("threadsNum", 't', "test problems dumension", false, 1,
     cmdline::range(1, 32));
  parser.add<double>("reliability", 'r', "reliability parameter for the method",
    false, 3.5, cmdline::range(1., 1000.));
  parser.add<double>("accuracy", 'e', "accuracy of the method",  false, 0.01);
  parser.add<int>("trialsLimit", 'l', "limit of trials for the method", false, 5000000);
  parser.parse_check(argc, argv);

  const unsigned nProblems = 100;
  ProblemsPool<gkls::GKLSFunction> pool;
  for(unsigned i = 0; i < nProblems; i++)
  {
    gkls::GKLSFunction* func = new gkls::GKLSFunction();
    func->SetFunctionClass(gkls::Simple, parser.get<int>("dimension"));
    func->SetType(gkls::TD);
    func->SetFunctionNumber(i + 1);
    pool.AddProblem(std::shared_ptr<gkls::GKLSFunction>(func));
  }

  std::cout << "Problems pool created\n";

  SolverParameters parameters(parser.get<double>("accuracy"),
      parser.get<double>("reliability"),
      parser.get<int>("threadsNum"),
      parser.get<int>("trialsLimit"), StopType::OptimumVicinity);
  parameters.logDeviations = false;

  GOSolver<gkls::GKLSFunction> solver;
  solver.SetParameters(parameters);
  solver.SetProblemsPool(pool);
  std::cout << "Solver started\n";
  auto start = std::chrono::system_clock::now();
  solver.Solve();
  auto end = std::chrono::system_clock::now();
  std::vector<Trial> optimumEstimations = solver.GetOptimumEstimations();
  std::vector<StatPoint> statistics = solver.GetStatistics();

  //for (size_t i = 0; i < statistics.size(); i++)
  //  std::cout << "(" << statistics[i].trial << ": " << statistics[i].meanDev << " " <<
  //   statistics[i].maxDev << ")\n";
  //for(size_t i = 0; i < optimumEstimations.size(); i++)
  //  std::cout << "Optimum value in problem #" << i + 1 << ": " << optimumEstimations[i].z << "\n";

  std::cout << "Number of trials: " << solver.GetTrialsNumber()
    << "\nNumber of iterations: " << solver.GetIterationsNumber() << "\n";

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Time elapsed: " << elapsed_seconds.count() << "s\n";

  return 0;
}
