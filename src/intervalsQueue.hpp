#ifndef OPTIMIZER_QUEUE_HPP
#define OPTIMIZER_QUEUE_HPP

#include "dataTypes.hpp"

class IntervalsQueue
{
protected:
	int MaxSize;
	int CurSize;
	Interval *pMem;
	int GetIndOfMinElem();
	void DeleteMinElem();
	void ReBuild(int Index);

public:
	// размер очереди должен быть равен 2^k - 1
	IntervalsQueue(int _MaxSize = 131071);
	~IntervalsQueue();
	
	int GetSize() const;
	bool IsEmpty() const;
	bool IsFull() const;

	void Push(const Interval &value);
	void PushWithPriority(const Interval &value);
	void DeleteInterval(const Interval &value);
	Interval Pop();

	void Clear();
};

#endif
