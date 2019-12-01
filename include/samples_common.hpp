#include <functional>
#include <cmath>
#include <chrono>
#include <iostream>

std::function<double()> buildComputeLoad(double delay)
{
  if (delay == 0)
    return [] {return 1;};

  double estimatedTime = 0;
  unsigned complexity = 0;
  unsigned delta = 500;

  auto computeKernel = [](unsigned iters)
  {
    double value = 0;
    for (unsigned i = 0; i < iters; i++)
    {
      double a1 = sin(value + i);
      double a2 = cos(value + i);
      value = a2*a2 + a1*a1;
    }
    return value + 1.;
  };

  double val = 0;
  do
  {
    complexity += delta;
    auto start = std::chrono::system_clock::now();
    for(int i = 0; i < 100; i++)
    {
      val = computeKernel(complexity);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    estimatedTime = elapsed_seconds.count() / 100;
  }
  while(estimatedTime * 1000. < delay);
  std::cout << "Estimated delay: " << estimatedTime*1000 << "\t complexity: " << complexity << "\t kernel value: " << val << '\n';

  return std::bind(computeKernel, complexity);
}
