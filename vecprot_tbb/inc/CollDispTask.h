#ifndef COLLDISPTASK
#define COLLDISPTASK

#include "tbb/task.h"
using namespace tbb;

#include "Rtypes.h"

class PropTask;

class CollDispTask : public task
{
private:
	Int_t fNumOfCollsToPop;				      // Number of collections to be poped by this task

public:
	CollDispTask (Int_t inNumOfCollsToPop);
	~CollDispTask ();

	task* execute ();

   PropTask& StartPropTasks (Int_t amountPriority, Int_t amountNormal);

};

#endif // COLLDISPTASK
