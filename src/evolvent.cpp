#include "evolvent.hpp"
#include <cassert>

Evolvent::Evolvent()
{
	mIsInitialized = false;
}

Evolvent::~Evolvent()
{
}

Evolvent::Evolvent(int dimension, int tightness, MapType type)
{
	assert(dimension > 1);
	assert(tightness > 2);

	mDimension = dimension;
	mTightness = tightness;
	mMapType = type;

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
