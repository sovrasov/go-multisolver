#include <json.hpp>
#include <cmdline.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <memory>
#include <omp.h>

#include <GKLS/GKLSConstrainedProblem.hpp>
#include <Grishagin/GrishaginConstrainedProblem.hpp>

#include "constrainedProblemsPool.hpp"
#include "goSolver.hpp"

int main(int argc, char** argv)
{
  cmdline::parser parser;
  parser.add<int>("dimension", 'd', "dimension of test problems class", false, 2,
     cmdline::range(1, 5));
  parser.add<int>("threadsNum", 'p', "number of threads in parallel method", false, 1,
     cmdline::range(1, 32));
  parser.add<int>("evolventTightness", 't', "", false, 12,
        cmdline::range(9, 16));
  parser.add<double>("reliability", 'r', "reliability parameter for the method",
    false, 4.5, cmdline::range(1., 1000.));
  parser.add<double>("accuracy", 'e', "accuracy of the method", false, 0.01);
  parser.add<int>("trialsLimit", 'l', "limit of trials for the method", false, 5000000);
  parser.add("saveStatistics", 's', "determines whether the method will "
    "save statistics of deviations");
  parser.add("verbose", 'v', "determines whether the method print information messages");
  parser.add<std::string>("statFile", 'f', "name of the file to write statistics",
    false, "statistics.csv");
  parser.add<std::string>("runMode", 'm', "", false, "multi", cmdline::oneof(
    std::string("multi"), std::string("synch"), std::string("asynch")));
  parser.add("hard", 'h', "determines type of GKLS functions");
  parser.add("mixedClass", 'c', "use mixed GKLS/Grishagin problems pool");
  parser.parse_check(argc, argv);

  GKLSClass problemsClass = parser.exist("hard") ? Hard : Simple;
  SolverParameters parameters(parser.get<double>("accuracy"),
    parser.get<double>("reliability"),
    parser.get<int>("threadsNum"),
    parser.get<int>("trialsLimit"), StopType::OptimumVicinity);
  parameters.verbose = parser.exist("verbose");
  parameters.evloventTightness = parser.get<int>("evolventTightness");
  std::vector<StatPoint> statistics;
  const unsigned nProblems = 100;

  auto start = std::chrono::system_clock::now();
  if (parser.get<std::string>("runMode") == std::string("multi"))
  {
    parameters.logDeviations = parser.exist("saveStatistics");
    GCGenProblemsPool<IConstrainedOptProblem> pool;

    for (unsigned i = 0; i < nProblems; i++)
    {
      if(i % 2 == 0 || !parser.exist("mixedClass"))
      {
        auto* problem = new TGKLSConstrainedProblem(cptInFeasibleDomain, 0.3, 0,
                                                    i + 1, parser.get<int>("dimension"), problemsClass, TD);
        pool.Add(std::shared_ptr<IConstrainedOptProblem>(problem));
      }
      else
      {
        auto* problem = new GrishaginConstrainedProblem(cptInFeasibleDomain, 0.3, 0, i + 1);
        pool.Add(std::shared_ptr<IConstrainedOptProblem>(problem));
      }
    }
    std::cout << "Problems pool created\n";

    GOSolver<GCGenProblemsPool<IConstrainedOptProblem>> solver;
    solver.SetParameters(parameters);
    solver.SetProblemsPool(pool);
    std::cout << "Solver started\n";
    solver.Solve();

    std::vector<Trial> optimumEstimations = solver.GetOptimumEstimations();
    statistics = solver.GetStatistics();

    std::cout << "Number of trials: " << solver.GetTrialsNumber()
      << "\nNumber of iterations: " << solver.GetIterationsNumber() << "\n";
    if(!statistics.empty())
      std::cout << "Problems solved: " << statistics.back().problems_solved << std::endl;
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
      GCGenProblemsPool<IConstrainedOptProblem> pool;
      IConstrainedOptProblem* problem;
      if(i % 2 == 0 || !parser.exist("mixedClass"))
      {
        problem = new TGKLSConstrainedProblem(cptInFeasibleDomain, 0.3, 0,
                                                    i + 1, parser.get<int>("dimension"), problemsClass, TD);
        pool.Add(std::shared_ptr<IConstrainedOptProblem>(problem));
      }
      else
      {
        problem = new GrishaginConstrainedProblem(cptInFeasibleDomain, 0.3, 0, i + 1);
        pool.Add(std::shared_ptr<IConstrainedOptProblem>(problem));
      }

      GOSolver<GCGenProblemsPool<IConstrainedOptProblem>> solver;
      solver.SetParameters(parameters);
      solver.SetProblemsPool(pool);
      solver.Solve();

      auto optPoint = problem->GetOptimumPoint();
      double currentDev = solver_internal::vectorsMaxDiff(optPoint.data(),
        solver.GetOptimumEstimations()[0].y, problem->GetDimension());

#pragma omp critical
      {
        StatPoint nextDeviation;
        nextDeviation.trial = solver.GetTrialsNumber() + lastStatistics.trial;
        if(statistics.size() != nProblems - 1)
          nextDeviation.maxDev = std::max(lastStatistics.maxDev, currentDev);
        else
          nextDeviation.maxDev = currentDev;
        if (currentDev < parameters.eps)
          nextDeviation.problems_solved = lastStatistics.problems_solved + 1;
        nextDeviation.meanDev = lastStatistics.meanDev + (-2. + currentDev) / 100.;

        statistics.push_back(nextDeviation);
        lastStatistics = nextDeviation;
        std::cout << "Problem # " << i + 1 <<
          " solved. Trials performed: " << solver.GetTrialsNumber() <<
          " Dev: " << nextDeviation.meanDev << "\n";
      }
    }
    std::cout << "Trials performed: " << lastStatistics.trial << "\n";
    std::cout << "Problems solved: " << lastStatistics.problems_solved << "\n";
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
