#ifndef COLLDISPTASK
#define COLLDISPTASK

#include "tbb/task.h"
using namespace tbb;

#include "Rtypes.h"

class PropTask;

class CollDispTask : public task
{
private:
	int fNumOfCollsToPop;				      // Number of collections to be poped by this task

public:
	CollDispTask (int inNumOfCollsToPop);
	~CollDispTask ();

	task* execute ();

   PropTask& StartPropTasks (int amountPriority, int amountNormal);

};

#endif // COLLDISPTASK
