#include "PropTask.h"

#include "GeantPropagator.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantVolumeBasket.h"
#include "GeantOutput.h"
#include "PhysicsProcess.h"
#include "Event_Hit.h"
#include "CollDispTask.h"

#include "TArrayI.h"
#include "TGeoManager.h"
#include "TGeoBranchArray.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TH1.h"

#include "tbb/task_scheduler_init.h"

#include <iostream>

PropTask::PropTask (Bool_t inPriority):
   fPriority (inPriority)
{ }

PropTask::~PropTask () { }

task* PropTask::execute ()
{
   GeantPropagator *propagator = GeantPropagator::Instance();
   tbb::task_scheduler_init init(propagator->fNthreads);
   PerThread::reference TBBperThread = propagator->fTBBthreadData.local();

   if (!fPriority) propagator->pnTasksTotal++;
   else propagator->ppTasksTotal++;

   Int_t npnTasks;
   Int_t nppTasks;
   /*Int_t iter3;
   Int_t iter4;*/
   if (propagator->fUseGraphics) {
	   if (!fPriority) {
         npnTasks = propagator->pnTasksRunning.fetch_and_increment();
         /*iter3 = propagator->niter3.fetch_and_increment();
         propagator->numOfnPropTasks->Fill(iter3, npnTasks+1);*/
	   } else {
         nppTasks = propagator->ppTasksRunning.fetch_and_increment();
         /*iter4 = propagator->niter4.fetch_and_increment();
         propagator->numOfpPropTasks->Fill(iter4, nppTasks+1);*/
      }
   }

   WorkloadManager *wm = WorkloadManager::Instance();
   Int_t *particles = TBBperThread.fPartInd->GetArray();
   Int_t *partnext  = TBBperThread.fPartNext->GetArray();
   Int_t *parttodo  = TBBperThread.fPartTodo->GetArray();
   Int_t *partcross = TBBperThread.fPartCross->GetArray();

   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) {
      nav = gGeoManager->AddNavigator();
      std::cerr << "[PropTask] Added navigator" << std::endl;
   }

   // Pop basket to propagate
   GeantBasket* basket;

   if (!fPriority) wm->tbb_feeder_queue.pop(basket);
   else wm->tbb_feeder_priority_queue.pop(basket);

   if (!basket) Fatal("PropTask", "No basket poped");

   Int_t ntotransport = basket->GetNtracks();
   if (!ntotransport) Fatal("PropTask", "Empty basket poped");

   if (propagator->fUseGraphics) {
	   if (fPriority) propagator->numOfTracksInPriorBasket->Fill (ntotransport);
	   else propagator->numOfTracksInBasket->Fill (ntotransport);
   }

   TBBperThread.fTracksPerBasket = ntotransport;

   GeantTrack **tracks = propagator->fTracks;
   TGeoBranchArray* path = tracks[basket->GetTracks()[0]]->path;
   TBBperThread.fVolume = path->GetCurrentNode()->GetVolume();
   GeantVolumeBasket* basket_sch = (GeantVolumeBasket*)TBBperThread.fVolume->GetField();

   memcpy(particles, basket->GetTracks(), ntotransport*sizeof(Int_t));

   GeantTrack *track = 0;
   Int_t ntotnext = 0;
   Int_t ntmp = 0;
   Int_t ntodo = 0;
   Int_t ncross = 0;
   Int_t cputime = 0;
   Int_t generation = 0;

   while (ntotransport)
   {
      generation++;
      if (propagator->fUsePhysics) propagator->PhysicsSelect(ntotransport, particles);
      basket_sch->ComputeTransportLength(ntotransport, particles);
      ntmp = ntotransport;
      ntodo = ntotransport;
      ntotnext = 0;
      Int_t *ptrParticles = particles;

      while (ntodo) {
         ntodo = 0;
         ncross = 0;
         // Propagate ALL ntmp tracks
         basket_sch->PropagateTracks(ntmp, ptrParticles, ntotnext, partnext, ntodo, parttodo, ncross, partcross);

         // Post-step actions by continuous processes for particles reaching boundaries
         if (propagator->fUsePhysics && ncross) {
            for (Int_t iproc=0; iproc < propagator->fNprocesses; iproc++) {
               if (propagator->Process(iproc)->IsType(PhysicsProcess::kDiscrete)) continue;
               Int_t nafter = 0;
               propagator->Process(iproc)->PostStep(TBBperThread.fVolume, ncross, partcross, nafter, NULL);
               basket_sch->ResetStep(ncross, partcross);
            }
         }
         ntmp = ntodo;
         ptrParticles = parttodo;
      }

      // Copy only tracks that survived boundaries (well we will have to think of
      // those too, like passing them to the next volume...)
      memcpy(particles, partnext, ntotnext*sizeof(Int_t));
      ntotransport = ntotnext;

      // Do post-step actions on remaining particles
      ntotnext = 0;
      // Loop all processes to group particles per process
      if (propagator->fUsePhysics && ntotransport) {
         // Apply continuous processes to all particles
         for (Int_t iproc=0; iproc < propagator->fNprocesses; iproc++) {
            if (propagator->Process(iproc)->IsType(PhysicsProcess::kDiscrete)) continue;
            ntodo = 0;
            propagator->Process(iproc)->PostStep(TBBperThread.fVolume, ntotransport, particles, ntodo, parttodo);
            // Do we have stopped particles ?
            if (ntodo<ntotransport) {
               memcpy(particles, parttodo, ntodo*sizeof(Int_t));
               ntotransport = ntodo;
            }
         }
         // Copy al tracks for which step was limited by a continuous process
         // to the next array
         for (Int_t itr=0; itr<ntotransport; itr++) {
            if (propagator->Process(tracks[particles[itr]]->process)->IsType(PhysicsProcess::kContinuous))
               partnext[ntotnext++] = particles[itr];
         }
         // Discrete processes only
         for (Int_t iproc=0; iproc < propagator->fNprocesses; iproc++) {
            // Make arrays of particles per process -> ntodo, parttodo
            if (propagator->Process(iproc)->IsType(PhysicsProcess::kContinuous)) continue;
            ntodo = 0;
            propagator->SelectTracksForProcess(iproc, ntotransport, particles, ntodo, parttodo);
            if (!ntodo) continue;
            if (TBBperThread.fPartTodo->GetSize()-ntodo<500) {
               TBBperThread.fPartTodo->Set(2*TBBperThread.fPartTodo->GetSize());
               parttodo  = TBBperThread.fPartTodo->GetArray();
            }
            // Do post step actions for particles suffering a given process.
            // Surviving particles are added to the next array
            propagator->Process(iproc)->PostStep(TBBperThread.fVolume, ntodo, parttodo, ntotnext, partnext);
            if (TBBperThread.fPartNext->GetSize()-ntotnext<500) {
               TBBperThread.fPartNext->Set(2*TBBperThread.fPartNext->GetSize());
               partnext  = TBBperThread.fPartNext->GetArray();
               TBBperThread.fPartInd->Set(2*TBBperThread.fPartInd->GetSize());
               particles = TBBperThread.fPartInd->GetArray();
            }
         }
         memcpy(particles, partnext, ntotnext*sizeof(Int_t));
         ntotransport = ntotnext;
      }
      // I/O: Dump current generation
//         Printf("   ### Generation %d:  %d tracks  cputime=%f", generation, ntotransport,cputime);
      if (propagator->fFillTree) {
         cputime = propagator->fTimer->CpuTime();
         propagator->fOutput->SetStamp(TBBperThread.fVolume->GetNumber(), generation, ntotransport, cputime);
         for (Int_t itrack=0; itrack<ntotransport;itrack++) {
            track = tracks[particles[itrack]];
            propagator->fOutput->SetTrack(itrack, track);
         }
         propagator->fOutTree->Fill();
      }
   }

   // Return empty basket
   basket->Clear();
   wm->tbb_feeder_empty_queue.push(basket);


//--------------------------------------------------------------------------------------------------------------------

   Int_t numOfFinishedEvents = TBBperThread.GetSizeOfFinishedEvents();
   Int_t numOfTracksInResultColl = TBBperThread.fCollection->GetNtracks();
   Int_t numOfCollPushed = 0;

   // Inject result collection if it is not empty
   if (numOfTracksInResultColl > 0) {
      propagator->InjectCollection (TBBperThread.fCollection);
      numOfCollPushed++;
   }

   // Check whether any priority events have finished
   if ((propagator->fPriorityRange[0] > -1) && numOfFinishedEvents) {
      Int_t first_not_transported = -1;
      for (Int_t evn=0; evn<propagator->fNtotal; evn++) {
         if (propagator->fEventsStatus[evn] == 0) {
            first_not_transported = evn;
            break;
         }
      }
      if (first_not_transported > propagator->fPriorityRange[1]) {
         Printf ("Stopping priority regime.");
         propagator->SetPriorityRange(-1, -1);
      } else if (first_not_transported== -1) {
         Printf ("first_not_transported = -1");
         Printf ("Stopping priority regime.");
         propagator->SetPriorityRange(-1, -1);
      }
   }

   // Check for finished events and inject new if needed
   if (numOfFinishedEvents > 0) {
      for (std::vector<Int_t>::iterator it = TBBperThread.fFinishedEvents.begin();
      it != TBBperThread.fFinishedEvents.end(); it++) {

         Int_t newEventInd = propagator->fNimportedEvents.fetch_and_increment();
         if (newEventInd < propagator->fNtotal) {
            numOfTracksInResultColl += propagator->ImportTracks(1, propagator->fNaverage, newEventInd, *it);
            numOfCollPushed++;
         } else break;

      }
   }

   TBBperThread.fFinishedEvents.clear();

   // Check whether we need to start dispatching
   Int_t nCollsToPop = 0;
   Int_t flagToStartDisp = kFALSE;

propagator->fPropTaskLock.Lock ();			// CRITICAL SECTION begin

   propagator->fCollsWaiting += numOfCollPushed;
   propagator->fTracksWaiting += numOfTracksInResultColl;

	// Non-priority regime
   if (propagator->fPriorityRange[0] == -1) {
	   if (propagator->fTracksWaiting >= propagator->fDispThr) {
		   nCollsToPop = propagator->fCollsWaiting;
		   propagator->fCollsWaiting = 0;
		   propagator->fTracksWaiting = 0;
		   flagToStartDisp = kTRUE;
	   }
   } else {
	      if ((propagator->fTracksWaiting >= propagator->fDispThrPriority) || (propagator->fGarbageCollMode)) {
		      nCollsToPop = propagator->fCollsWaiting;
		      propagator->fCollsWaiting = 0;
		      propagator->fTracksWaiting = 0;
		      flagToStartDisp = kTRUE;
	      }
   }

propagator->fPropTaskLock.UnLock ();		// CRITICAL SECTION end

   // Start dispatching task
	if (flagToStartDisp) {
		empty_task& cont = *new(task::allocate_continuation()) empty_task();
		cont.set_ref_count(1);
		CollDispTask& dispTask = *new(cont.allocate_child()) CollDispTask(nCollsToPop);

      if (propagator->fUseGraphics) {
	      if (!fPriority) {
            npnTasks = propagator->pnTasksRunning.fetch_and_decrement();
            /*iter3 = propagator->niter3.fetch_and_increment();
            propagator->numOfnPropTasks->Fill(iter3, npnTasks-1);*/
	      } else {
            nppTasks = propagator->ppTasksRunning.fetch_and_decrement();
            /*iter4 = propagator->niter4.fetch_and_increment();
            propagator->numOfpPropTasks->Fill(iter4, nppTasks-1);*/
         }
      }

		return &dispTask;
	}

   if (propagator->fUseGraphics) {
      if (!fPriority) {
         npnTasks = propagator->pnTasksRunning.fetch_and_decrement();
         /*iter3 = propagator->niter3.fetch_and_increment();
         propagator->numOfnPropTasks->Fill(iter3, npnTasks-1);*/
      } else {
         nppTasks = propagator->ppTasksRunning.fetch_and_decrement();
         /*iter4 = propagator->niter4.fetch_and_increment();
         propagator->numOfpPropTasks->Fill(iter4, nppTasks-1);*/
      }
   }

   return NULL;
}

