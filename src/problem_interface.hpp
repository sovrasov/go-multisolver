#pragma once

class IGOPRoblem
{
public:
	virtual double Calculate(const double* y) const = 0;
	virtual int GetOptimumCoordinates(double* y) const = 0;
	virtual void GetDomainBounds(double* lb, double* ub) const = 0;
	virtual unsigned GetDimension() const = 0;
	virtual double GetOptimalValue() const = 0;
};
