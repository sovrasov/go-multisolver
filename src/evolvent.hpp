#ifndef EVOLVENT_HPP
#define EVOLVENT_HPP

#include "Map.hpp"

#define MAX_PREIMAGES 32

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
	Evolvent(int dimension, int tightness, MapType type = Simple);
	~Evolvent();

	void GetImage(double x, double y[]);
	int GetAllPreimages(double* p, double xp[]);
};


#endif
