#include "CollDispTask.h"

#include "GeantPropagator.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "PropTask.h"
#include "GeantMainScheduler.h"

#include "TH1.h"

#include "tbb/task_scheduler_init.h"

CollDispTask::CollDispTask (Int_t inNumOfCollsToPop):
	fNumOfCollsToPop (inNumOfCollsToPop)
{ }

CollDispTask::~CollDispTask () { }

task* CollDispTask::execute ()
{
   GeantPropagator *propagator = GeantPropagator::Instance();
   PerThread::reference TBBperThread = propagator->fTBBthreadData.local();
   WorkloadManager *wm = WorkloadManager::Instance();

   propagator->dTasksRunning++;

   tbb::task_scheduler_init init(propagator->fNthreads+1);

   // Number of baskets in the queue to transport
   Int_t ntotransport = 0;

   UInt_t npop = 0;
   Int_t ninjectedNormal = 0;
   Int_t ninjectedPriority = 0;

   GeantTrackCollection *collector;
   TObject **carray = new TObject*[500];

   // Pop certain number of collections recieved when constructing task object
   for (Int_t i=0; i<fNumOfCollsToPop; i++) {
      wm->tbb_collector_queue.pop(collector);
      carray[i] = (TObject*)collector;
      npop++;
   }

   propagator->numOfCollsPerTask->Fill(npop);

   // Process popped collections and flush their tracks
   for (UInt_t icoll=0; icoll<npop; icoll++) {
      collector = (GeantTrackCollection*)carray[icoll];

      propagator->numOfTracksInColl->Fill(collector->GetNtracks());

      collector->FlushTracks(wm->fMainSch, &ninjectedNormal, &ninjectedPriority);
      wm->tbb_collector_empty_queue.push(collector);
   }

   // Flush priority baskets if in priority mode
   if (propagator->fPrioritize) {
      ninjectedPriority += wm->fMainSch->FlushPriorityBaskets();
//      Printf ("ninjectedPriority=%d", ninjectedPriority);
   }

   ntotransport = wm->tbb_feeder_queue.size();

   // Low-populated feeder_queue - do something!
   if (ntotransport < propagator->fMinFeeder) {

      ninjectedNormal += wm->fMainSch->FlushNormalBaskets();

      if (!propagator->fPrioritize) {

         // Find first not transported
         propagator->fDispTaskLock.Lock();			// CRITICAL SECTION begin

         Int_t first_not_transported = -1;
         for (Int_t evn=0; evn<propagator->fNtotal; evn++) {
            if (propagator->fEventsStatus[evn] == 0) {
               first_not_transported = evn;
               break;
            }
         }

         // Prioritize
         if (first_not_transported > -1) {
            propagator->SetPriorityRange(first_not_transported,
               TMath::Min(first_not_transported+propagator->fNevToPrioritize-1, propagator->fNtotal-1));
            propagator->fPrioritize = kTRUE;
            Printf ("Prioritizing events %d to %d", first_not_transported,
               TMath::Min(first_not_transported+propagator->fNevToPrioritize-1, propagator->fNtotal-1));

         } else Printf ("Trying to find an event to prioritize when all the events are transported?");

         propagator->fDispTaskLock.UnLock();		// CRITICAL SECTION end

      }

   }

/*------------------------------------------------------------------------------------------------------*/
// Just histograms filling
///*
   Int_t localNiter = propagator->niter.fetch_and_increment();
   Double_t nperbasket = 0;
   for (PerThread::iterator it=propagator->fTBBthreadData.begin(); it!= propagator->fTBBthreadData.end(); ++it)
      nperbasket += it->fTracksPerBasket;
   nperbasket /= propagator->fNthreads+1;
   if (propagator->fUseGraphics) {
      if ((localNiter%10==0) && (localNiter/10<200000)) {
         propagator->hnb->Fill(localNiter/10, ntotransport);
         propagator->htracks->Fill(localNiter/10, nperbasket);
         propagator->numOfTracksTransportedInIter->Fill(localNiter/10, propagator->fNtransported);
      }
   }
//*/
/*------------------------------------------------------------------------------------------------------*/

	if (ninjectedNormal+ninjectedPriority) {
		PropTask& lastTask = StartPropTasks (ninjectedPriority, ninjectedNormal);
      propagator->dTasksRunning--;
		return &lastTask;
	}

   propagator->dTasksRunning--;
   return NULL;
}

// total amount must be more then 0
PropTask& CollDispTask::StartPropTasks (Int_t amountPriority, Int_t amountNormal)
{
	task_list tlist;
	empty_task& cont = *new(allocate_continuation()) empty_task();
	cont.set_ref_count (amountPriority+amountNormal);
	Bool_t flag;

	if (amountPriority) {

		for (Int_t i = 0; i < amountPriority - 1; i++)
			tlist.push_back (*new(cont.allocate_child()) PropTask(kTRUE));
		for (Int_t i = 0; i < amountNormal; i++)
			tlist.push_back (*new(cont.allocate_child()) PropTask(kFALSE));
		flag = kTRUE;

	} else if (amountNormal) {

		for (Int_t i = 0; i < amountPriority; i++)
			tlist.push_back (*new(cont.allocate_child()) PropTask(kTRUE));
		for (Int_t i = 0; i < amountNormal-1; i++)
			tlist.push_back (*new(cont.allocate_child()) PropTask(kFALSE));
		flag = kFALSE;
	}

	PropTask& lastTask = *new(cont.allocate_child()) PropTask(flag);

	if (!tlist.empty()) spawn(tlist);
	return lastTask;
}

