#include "Grishagin/grishagin_function.hpp"
#include "GKLS/gkls_function.hpp"
#include "problemsPool.hpp"
#include "goSolver.hpp"

#include <iostream>
#include <fstream>
#include <chrono>
#include <memory>
#include <cmdline.h>
#include <omp.h>

int main(int argc, char** argv)
{
  cmdline::parser parser;
  parser.add<int>("dimension", 'd', "dimension of test problems class", false, 2,
     cmdline::range(1, 5));
  parser.add<int>("threadsNum", 't', "number of threads in parallel method", false, 1,
     cmdline::range(1, 32));
  parser.add<double>("reliability", 'r', "reliability parameter for the method",
    false, 4.5, cmdline::range(1., 1000.));
  parser.add<double>("accuracy", 'e', "accuracy of the method", false, 0.01);
  parser.add<int>("trialsLimit", 'l', "limit of trials for the method", false, 5000000);
  parser.add("saveStatistics", 's', "determines whether the method will "
    "save statistics of deviations");
  parser.add<std::string>("statFile", 'f', "name of the file to write statistics",
    false, "statistics.csv");
  parser.add<std::string>("runMode", 'm', "", false, "multi", cmdline::oneof(
    std::string("multi"), std::string("synch"), std::string("asynch")));
  parser.add("hard", 'h', "determines type of GKLS functions");
  parser.parse_check(argc, argv);

  gkls::GKLSClass problemsClass = parser.exist("hard") ? gkls::Hard : gkls::Simple;
  SolverParameters parameters(parser.get<double>("accuracy"),
    parser.get<double>("reliability"),
    parser.get<int>("threadsNum"),
    parser.get<int>("trialsLimit"), StopType::OptimumVicinity);
  std::vector<StatPoint> statistics;
  const unsigned nProblems = 100;

  auto start = std::chrono::system_clock::now();
  if (parser.get<std::string>("runMode") == std::string("multi"))
  {
    parameters.logDeviations = parser.exist("saveStatistics");
    ProblemsPool<gkls::GKLSFunction> pool;
    for (unsigned i = 0; i < nProblems; i++)
    {
      gkls::GKLSFunction* func = new gkls::GKLSFunction();
      func->SetFunctionClass(problemsClass, parser.get<int>("dimension"));
      func->SetType(gkls::TD);
      func->SetFunctionNumber(i + 1);
      pool.Add(std::shared_ptr<gkls::GKLSFunction>(func));
    }

    std::cout << "Problems pool created\n";

    GOSolver<gkls::GKLSFunction> solver;
    solver.SetParameters(parameters);
    solver.SetProblemsPool(pool);
    std::cout << "Solver started\n";
    solver.Solve();
    std::vector<Trial> optimumEstimations = solver.GetOptimumEstimations();
    statistics = solver.GetStatistics();

    std::cout << "Number of trials: " << solver.GetTrialsNumber()
      << "\nNumber of iterations: " << solver.GetIterationsNumber() << "\n";
    //  for(size_t i = 0; i < optimumEstimations.size(); i++)
    //    std::cout << "Optimum value in problem #" << i + 1 << ": " << optimumEstimations[i].z << "\n";
  }
  else
  {
    parameters.numThreads = 1;
    StatPoint lastStatistics(0, 2., 2.);

#if defined _MSC_VER

    if(parser.get<std::string>("runMode") == std::string("asynch"))
    {
      std::cerr << "asynch mode is disabled on current platform\n";
      return 0;
    }

#pragma omp parallel for num_threads(parser.get<int>("threadsNum")), schedule(static, 1)
#else

    if(parser.get<std::string>("runMode") == std::string("synch"))
      omp_set_schedule(omp_sched_static, 1);
    else if(parser.get<std::string>("runMode") == std::string("asynch"))
      omp_set_schedule(omp_sched_dynamic, 1);

#pragma omp parallel for num_threads(parser.get<int>("threadsNum")), schedule(runtime)
#endif
    for (int i = 0; i < (int)nProblems; i++)
    {
      ProblemsPool<gkls::GKLSFunction> pool;
      gkls::GKLSFunction* func = new gkls::GKLSFunction();
      func->SetFunctionClass(problemsClass, parser.get<int>("dimension"));
      func->SetType(gkls::TD);
      func->SetFunctionNumber(i + 1);
      pool.Add(std::shared_ptr<gkls::GKLSFunction>(func));

      GOSolver<gkls::GKLSFunction> solver;
      solver.SetParameters(parameters);
      solver.SetProblemsPool(pool);
      solver.Solve();

      double optPoint[solverMaxDim];
      func->GetOptimumCoordinates(optPoint);
      double currentDev = solver_internal::vectorsMaxDiff(optPoint,
        solver.GetOptimumEstimations()[0].y, func->GetDimension());

#pragma omp critical
      {
        StatPoint nextDeviation;
        nextDeviation.trial = solver.GetTrialsNumber() + lastStatistics.trial;
        if(statistics.size() != nProblems - 1)
          nextDeviation.maxDev = std::max(lastStatistics.maxDev, currentDev);
        else
          nextDeviation.maxDev = currentDev;
        nextDeviation.meanDev = lastStatistics.meanDev + (-2. + currentDev) / 100.;

        statistics.push_back(nextDeviation);
        lastStatistics = nextDeviation;
        std::cout << "Problem # " << i + 1 <<
          " solved. Trials performed: " << solver.GetTrialsNumber() <<
          " Dev: " << nextDeviation.meanDev << "\n";
      }
    }
    std::cout << "Trials performed: " << lastStatistics.trial << "\n";
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Time elapsed: " << elapsed_seconds.count() << "s\n";

  if(parser.exist("saveStatistics"))
  {
    std::ofstream fout;
    fout.open(parser.get<std::string>("statFile"), std::ios_base::out);

    for (size_t i = 0; i < statistics.size(); i++)
      fout << statistics[i].trial << ", " << statistics[i].meanDev << ", " <<
       statistics[i].maxDev << "\n";

   fout.close();
  }


  return 0;
}
