#pragma once

#include <cstring>

const unsigned solverMaxDim = 5;
const unsigned solverMaxFunctionsNum = 5;

struct Trial
{
  double x;
  double y[solverMaxDim];
  int v;
  double z[solverMaxFunctionsNum];
  Trial() : v(0) {}
  Trial(double _x) : x(_x), v(0) {}

  double GetZ() const
  {
    return z[v];
  }
};

struct Point
{
  double x;
  int v;
  double z[solverMaxFunctionsNum];
  Point() : v(0) {}
  Point(double _x) : x(_x), v(0) {}
  Point(const Trial& t) : x(t.x), v(t.v)
  {
    std::memcpy(z, t.z, solverMaxFunctionsNum*sizeof(double));
  }

  double GetZ() const
  {
    return z[v];
  }
};

struct Interval
{
  Point xl;
  Point xr;
  double R;
  double delta;
  unsigned problemIdx;
  Interval() {}
  Interval(const Point& _xl, const Point& _xr) : xl(_xl), xr(_xr) {}
};

struct CompareIntervals
{
  bool operator() (const Interval* i1, const Interval* i2) const
  {
    return i1->xl.x < i2->xl.x;
  }
};

class CompareByR
{
public:
  bool operator() (const Interval* i1, const Interval* i2)
  {
    return i1->R < i2->R;
  }
};

struct StatPoint
{
  unsigned trial;
  double maxDev;
  double meanDev;
  unsigned problems_solved;

  StatPoint() {}
  StatPoint(unsigned _trial, double _maxDev = 0., double _meanDev = 0.) :
    trial(_trial), maxDev(_maxDev), meanDev(_meanDev), problems_solved(0) {}
};
