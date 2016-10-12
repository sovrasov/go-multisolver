#include "evolvent.hpp"
#include <cassert>
#include <cmath>

Evolvent::Evolvent()
{
  mIsInitialized = false;
}

Evolvent::~Evolvent()
{
}

Evolvent::Evolvent(int dimension, int tightness, double* lb, double* ub, MapType type)
{
  assert(dimension > 1);
  assert(tightness > 2);

  mDimension = dimension;
  mTightness = tightness;
  mMapType = type;

  mShiftScalars.resize(mDimension);
  mRho = 0.;
  for (int i = 0; i < mDimension; i++)
	{
		mRho = fmax(mRho, ub[i] - lb[i]);
		mShiftScalars[i] = 0.5*(lb[i] + ub[i]);
	}

  switch (mMapType)
  {
  case Simple:
    mMapKey = 1;
    break;
  case Linear:
    mMapKey = 2;
    break;
  case Noninjective:
    mMapKey = 3;
    break;
  }

  mIsInitialized = true;
}

void Evolvent::GetImage(double x, double y[])
{
  mapd(x, mTightness, y, mDimension, mMapKey);

  for (int i = 0; i < mDimension; i++)
    y[i] = mRho*y[i] + mShiftScalars[i];
}

int Evolvent::GetAllPreimages(double * p, double xp[])
{
  int preimNumber = 1;
  if(mMapType == Noninjective)
    invmad(mTightness, xp, MAX_PREIMAGES, &preimNumber, p, mDimension, 4);
  else
    xyd(xp, mTightness, p, mDimension);

  return preimNumber;
}
