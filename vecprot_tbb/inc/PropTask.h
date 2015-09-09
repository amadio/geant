#ifndef PROPTASK
#define PROPTASK

#include "tbb/task.h"
using namespace tbb;

#include "Rtypes.h"


class PropTask : public task
{
private:
	bool fPriority;             // True if this task propagates a priority basket

public:
	PropTask (bool inPriority);
	~PropTask ();

	task* execute ();

};

#endif // PROPTASK
