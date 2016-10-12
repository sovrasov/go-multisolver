#ifndef OPTIMIZER_QUEUE_HPP
#define OPTIMIZER_QUEUE_HPP

#include "dataTypes.hpp"
#include <vector>

class IntervalsQueue
{
protected:
	int MaxSize;
	int CurSize;
	std::vector<Interval*> pMem;
	int GetIndOfMinElem();
	void DeleteMinElem();
	void ReBuild(int Index);

public:
	// размер очереди должен быть равен 2^k - 1
	IntervalsQueue(int _MaxSize = 262143);
	~IntervalsQueue();

	int GetSize() const;
	bool IsEmpty() const;
	bool IsFull() const;

	void Push(Interval* value);
	void PushWithPriority(Interval* value);
	Interval* Pop();

	void Clear();
};

#endif
