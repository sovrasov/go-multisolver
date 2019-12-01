#include <cmdline.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <omp.h>

#include "mco_problems_adapter.hpp"
#include "goSolver.hpp"
#include "samples_common.hpp"

bool isVectorLess(const double* v1, const double* v2, int dim, double filterEps = 0)
{
  for (int i = 0; i < dim; i++)
    if (v1[i] - filterEps >= v2[i])
      return false;

  return true;
}

std::vector<Trial> GetWeakOptimalPoints(const std::vector<Trial>& source_points, int num_objs, int num_constrs)
{
  std::vector<Trial> optTrials;
  const size_t data_size = source_points.size();

  for(size_t i = 0; i < data_size; i++)
  {
    if(source_points[i].v == num_constrs)
    {
      bool isWeakOptimal = true;
      for(size_t j = 0; j < data_size; j++)
      {
        if(i != j && source_points[j].v == num_constrs)
        {
          if(isVectorLess(source_points[j].z, source_points[i].z,
                  num_objs, 0.0))
            isWeakOptimal = false;
        }
      }
      if(isWeakOptimal)
        optTrials.push_back(source_points[i]);
    }
  }
  return optTrials;
}

int main(int argc, char** argv)
{
  cmdline::parser parser;
  parser.add<int>("threadsNum", 'p', "number of threads in parallel method", false, 1,
     cmdline::range(1, 32));
  parser.add<int>("evolventTightness", 't', "", false, 12,
        cmdline::range(9, 16));
  parser.add<double>("reserves", 0, "eps-reservation parameter", false, 0);
  parser.add<double>("reliability", 'r', "reliability parameter for the method",
    false, 4.5, cmdline::range(1., 1000.));
  parser.add<double>("accuracy", 'e', "accuracy of the method", false, 0.01);
  parser.add<int>("trialsLimit", 'l', "limit of trials for the method", false, 5000000);
  parser.add<unsigned>("problemsNum", 'n', "number of pareto points", false, 100);
  parser.add("saveStatistics", 's', "determines whether the method will "
    "save statistics of deviations");
  parser.add("verbose", 'v', "determines whether the method print information messages");
  parser.add<std::string>("statFile", 'f', "name of the file to write statistics",
    false, "statistics.csv");
  parser.add<std::string>("save_points", 0, "name of the file to save solution",
    false, "");
  parser.add<double>("delay", 0, "Delay in objective functions (ms)", false, 0);
  parser.add<std::string>("runMode", 'm', "", false, "multi", cmdline::oneof(
    std::string("multi"), std::string("synch"), std::string("asynch")));
  parser.parse_check(argc, argv);

  auto f1 = [](const double* x) {return 4*(x[0]*x[0] + x[1]*x[1]);};
  auto f2 = [](const double* x) {return pow(x[0] - 5, 2) + pow(x[1] - 5, 2);};
  auto g1 = [](const double* x) {return pow(x[0] - 5, 2) + x[1] * x[1] - 25;};
  auto g2 = [](const double* x) {return -pow(x[0] - 8, 2) - pow(x[1] + 3, 2) + 7.7;};

  auto computeLoad = buildComputeLoad(parser.get<double>("delay"));

  unsigned num_pareto_points = parser.get<unsigned>("problemsNum");
  SolverParameters parameters(parser.get<double>("accuracy"),
    parser.get<double>("reliability"),
    parser.get<int>("threadsNum"),
    parser.get<int>("trialsLimit"), StopType::Accuracy);
  parameters.verbose = parser.exist("verbose");
  parameters.evloventTightness = parser.get<int>("evolventTightness");
  parameters.epsR = parser.get<double>("reserves");

  std::vector<Trial> optimumEstimations;
  std::cout << "Solver started\n";
  auto start = std::chrono::system_clock::now();

  if (parser.get<std::string>("runMode") == std::string("multi"))
  {
    MultiObjectiveProblemAdapter problem({f1, f2}, {g1, g2}, {0, 0}, {5, 3}, num_pareto_points);
    problem.SetComputeLoad(computeLoad);
    GOSolver<MultiObjectiveProblemAdapter> solver;
    solver.SetParameters(parameters);
    solver.SetProblemsPool(problem);
    solver.Solve();
    optimumEstimations = solver.GetOptimumEstimations();
    std::cout << "Number of trials: " << solver.GetTrialsNumber()
      << "\nNumber of iterations: " << solver.GetIterationsNumber() << "\n";
  }
  else
  {
    parameters.numThreads = 1;
    parameters.trialsLimit /= num_pareto_points;
    double h = 1. / (num_pareto_points - 1);
    unsigned total_trials = 0;
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
    for (int i = 0; i < (int)num_pareto_points; i++)
    {
      MultiObjectiveProblemAdapter problem({f1, f2}, {g1, g2}, {0, 0}, {5, 3}, 1);
      problem.SetLambdas({i*h});
      problem.SetComputeLoad(computeLoad);
      GOSolver<MultiObjectiveProblemAdapter> solver;
      solver.SetParameters(parameters);
      solver.SetProblemsPool(problem);
      solver.Solve();
#pragma omp critical
      {
        optimumEstimations.push_back(solver.GetOptimumEstimations()[0]);
        total_trials += solver.GetTrialsNumber();
      }
    }
    std::cout << "Number of trials: " << total_trials << "\n";
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Time elapsed: " << elapsed_seconds.count() << "s\n";
  for (auto& paretoPoint : optimumEstimations)
  {
    paretoPoint.z[0] = f1(paretoPoint.y);
    paretoPoint.z[1] = f2(paretoPoint.y);
  }
  optimumEstimations = GetWeakOptimalPoints(optimumEstimations, 2, 2);

  auto save_points_path = parser.get<std::string>("save_points");
  std::ofstream points_file;
  if (!save_points_path.empty())
  {
    points_file.open(save_points_path);
    points_file << 2 << "\n" << 2 << "\n";
  }

  for (const auto& paretoPoint : optimumEstimations)
  {
    if (paretoPoint.v == 2)
    {
      //std::cout << "Pareto point: " << f1_val << ", " << f2_val << "\n";
      if (points_file.is_open())
        points_file << paretoPoint.y[0] << ", " << paretoPoint.y[1] << ", "
                    << paretoPoint.z[0] << ", " << paretoPoint.z[1] << ";\n";
    }
    else
      std::cout << "Failed to find a feasible point!" << '\n';
  }
  return 0;
}
