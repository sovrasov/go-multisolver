#ifndef EVOLVENT_HPP
#define EVOLVENT_HPP

#include "Map.hpp"
#include <vector>

#define MAX_PREIMAGES 32

class Evolvent
{
protected:
	int mDimension;
	int mTightness;
	
	double mRho;
	std::vector<double> mShiftScalars;

	bool mIsInitialized;
private:
	MapType mMapType;
	int mMapKey;

public:
	Evolvent();
	Evolvent(int dimension, int tightness, double* lb, double*ub, MapType type = Simple);
	Evolvent(Evolvent& source);
	~Evolvent();

	void GetImage(double x, double y[]);
	int GetAllPreimages(double* p, double xp[]);
};


#endif
