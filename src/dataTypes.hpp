#ifndef DATA_TYPES_H
#define DATA_TYPES_H

const unsigned solverMaxDim = 5;

struct Trial
{
  double x;
  double y[solverMaxDim];
  double z;
};

struct Interval
{
  Trial xl;
  Trial xr;
  double R;
  double delta;
  unsigned problemIdx;
};

#endif
