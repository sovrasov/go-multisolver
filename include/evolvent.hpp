#pragma once

#define MAX_PREIMAGES 32

enum class MapType {
  Simple = 1, Linear = 2, Noninjective = 3
};

class Evolvent
{
protected:
  int mDimension;
  int mTightness;
  bool mIsInitialized;

private:
  MapType mMapType;
  int mMapKey;

public:
  Evolvent();
  Evolvent(int dimension, int tightness, MapType type = MapType::Simple);
  ~Evolvent();

  void GetImage(double x, double y[]) const;
  void GetImage(double x, double y[], const double lb[], const double ub[]) const;
  int GetAllPreimages(double* p, double xp[]);
};
