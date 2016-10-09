#include "Grishagin/grishagin_function.hpp"
#include "GKLS/gkls_function.hpp"
#include "problemsPool.hpp"

#include <iostream>

int main(int argc, const char** argv)
{
  ProblemsPool<gkls::GKLSFunction> pool;
  gkls::GKLSFunction** functions = new gkls::GKLSFunction*[100];
  for(int i = 0; i < 100; i++)
  {
    gkls::GKLSFunction* func = new gkls::GKLSFunction();
    func->SetFunctionClass(gkls::Hard, 3);
    func->SetType(gkls::TD);
    func->SetFunctionNumber(1);
    pool.AddProblem(func);
  }

  std::cout << "Problems pool created\n";

  for(int i = 0; i < 100; i++)
    delete functions[i];
  delete[] functions;

  return 0;
}
