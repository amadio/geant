#ifndef PROPTASK
#define PROPTASK

#include "tbb/task.h"
using namespace tbb;

#include "Rtypes.h"



class PropTask : public task
{
private:
	Bool_t fPriority;             // True if this task propagates a priority basket

public:
	PropTask (Bool_t inPriority);
	~PropTask ();

	task* execute ();

};

#endif // PROPTASK
