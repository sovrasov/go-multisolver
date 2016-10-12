#include "intervalsQueue.hpp"

#include <cmath>
#include <algorithm>

IntervalsQueue::IntervalsQueue(int _MaxSize)
{
	MaxSize = _MaxSize;
	CurSize = 0;
  pMem.resize(MaxSize);// = new Interval*[MaxSize];
}

// ------------------------------------------------------------------------------------------------
IntervalsQueue::~IntervalsQueue()
{
	//delete[] pMem;
}
// ------------------------------------------------------------------------------------------------
bool IntervalsQueue::IsEmpty() const
{
	return CurSize == 0;
}
// ------------------------------------------------------------------------------------------------
int IntervalsQueue::GetSize() const
{
	return CurSize;
}

// ------------------------------------------------------------------------------------------------
bool IntervalsQueue::IsFull() const
{
	return CurSize == MaxSize;
}

// ------------------------------------------------------------------------------------------------
void IntervalsQueue::Push(Interval* value)
{
	if (IsFull())
	{
		int MinInd = GetIndOfMinElem();
		if (value->R > pMem[MinInd]->R)
			DeleteMinElem();
		else
			return;
	}
	CurSize++;
	pMem[CurSize - 1] = value;
	if (CurSize > 1)
		ReBuild(CurSize - 1);
}

// ------------------------------------------------------------------------------------------------
void IntervalsQueue::PushWithPriority(Interval* value)
{
	if (IsEmpty())
	{
		CurSize++;
		pMem[CurSize - 1] = value;
	}
	else
	{
		int MinInd = GetIndOfMinElem();

		if (value->R >= pMem[MinInd]->R)
		{
			if (IsFull())
				DeleteMinElem();
			CurSize++;
			pMem[CurSize - 1] = value;
			if (CurSize > 1)
				ReBuild(CurSize - 1);
		}
	}
}

// ------------------------------------------------------------------------------------------------
Interval* IntervalsQueue::Pop()
{
	Interval* tmp = pMem[0];
	pMem[0] = pMem[CurSize - 1];
	CurSize--;
	if (CurSize > 1)
		ReBuild(0);

	return tmp;
}

// ------------------------------------------------------------------------------------------------
int IntervalsQueue::GetIndOfMinElem()
{
	int i, StartIndex;
	double min = HUGE_VAL;
	int MinIndex = -1;

	if (CurSize % 2)
		StartIndex = (CurSize - 1) / 2;
	else
		StartIndex = (CurSize - 1) / 2 + 1;

	for (i = StartIndex; i < CurSize; i++)
		if (min > pMem[i]->R)
		{
			MinIndex = i;
			min = pMem[i]->R;
		}

	return MinIndex;
}

// ------------------------------------------------------------------------------------------------
void IntervalsQueue::DeleteMinElem()
{
	int MinInd = GetIndOfMinElem();
	pMem[MinInd] = pMem[CurSize - 1];
	CurSize--;
	if (CurSize > 1)
		ReBuild(MinInd);
}

// ------------------------------------------------------------------------------------------------
void IntervalsQueue::ReBuild(int Index)
{
	int i, j, k;
	if (Index == 0)
	{
		i = Index;
		j = 2 * i + 1;
		k = 2 * i + 2;
		if (k < CurSize)
			if (pMem[j]->R < pMem[k]->R)
				j = k;
		while (true)
		{
			if (pMem[i]->R >= pMem[j]->R)
				break;
			std::swap(pMem[i], pMem[j]);

			i = j;
			j = 2 * i + 1;
			k = 2 * i + 2;
			if (j > CurSize - 1)
				break;
			if (k < CurSize)
				if (pMem[j]->R < pMem[k]->R)
					j = k;
		}
	}
	else
	{
		i = Index;
		j = (i - 1) / 2;
		while ((i > 0) && (pMem[j]->R <= pMem[i]->R))
		{
			std::swap(pMem[i], pMem[j]);
			i = j;
			j = (i - 1) / 2;
		}
	}
}

void IntervalsQueue::Clear()
{
	CurSize = 0;
}
// - end of file ----------------------------------------------------------------------------------
