#include <json.hpp>
#include <cmdline.h>

#include <iostream>
#include <fstream>
#include <chrono>

#include "mco_problems_adapter.hpp"
#include "goSolver.hpp"

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
  parser.add<std::string>("runMode", 'm', "", false, "multi", cmdline::oneof(
    std::string("multi"), std::string("synch"), std::string("asynch")));
  parser.parse_check(argc, argv);

  auto f1 = [](const double* x) {return 4*(x[0]*x[0] + x[1]*x[1]);};
  auto f2 = [](const double* x) {return pow(x[0] - 5, 2) + pow(x[1] - 5, 2);};
  auto g1 = [](const double* x) {return pow(x[0] - 5, 2) + x[1] * x[1] - 25;};
  auto g2 = [](const double* x) {return -pow(x[0] - 8, 2) - pow(x[1] + 3, 2) + 7.7;};

  MultiObjectiveProblemAdapter problem({f1, f2}, {g1, g2}, {0, 0}, {5, 3}, 100);

  GOSolver<MultiObjectiveProblemAdapter> solver;
  SolverParameters parameters(parser.get<double>("accuracy"),
    parser.get<double>("reliability"),
    parser.get<int>("threadsNum"),
    parser.get<int>("trialsLimit"), StopType::Accuracy);
  parameters.verbose = parser.exist("verbose");
  parameters.evloventTightness = parser.get<int>("evolventTightness");
  parameters.epsR = parser.get<double>("reserves");

  solver.SetParameters(parameters);
  solver.SetProblemsPool(problem);
  std::cout << "Solver started\n";
  auto start = std::chrono::system_clock::now();
  solver.Solve();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Time elapsed: " << elapsed_seconds.count() << "s\n";
  std::cout << "Number of trials: " << solver.GetTrialsNumber()
    << "\nNumber of iterations: " << solver.GetIterationsNumber() << "\n";
  std::vector<Trial> optimumEstimations = solver.GetOptimumEstimations();

  auto save_points_path = parser.get<std::string>("save_points");
  std::ofstream points_file;
  if (!save_points_path.empty())
  {
    points_file.open(save_points_path);
    points_file << 2 << "\n" << 2 << "\n";
  }

  for (const auto& paretoPoint : optimumEstimations)
  {
    double f1_val = f1(paretoPoint.y);
    double f2_val = f2(paretoPoint.y);
    if (paretoPoint.v == 2)
    {
      //std::cout << "Pareto point: " << f1_val << ", " << f2_val << "\n";
      if (points_file.is_open())
        points_file << paretoPoint.y[0] << ", " << paretoPoint.y[1] << ", "
                    << f1_val << ", " << f2_val << ";\n";
    }
    else
      std::cout << "Failed to find a feasible point!" << '\n';
  }
  return 0;
}
