#ifndef GEANT_TASKMANAGER
#define GEANT_TASKMANAGER

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#include <vector>
#include "sync_objects.h"


typedef void (*ExecuteTaskFunc_t)(void *data);

//______________________________________________________________________________
//   GeantTaskManager - A general concurrent task manager steered by data
//                      availability. Use: RegisterFunction() method to register 
//                      an arbitrary number of functions of type: 
//                         void Function(void *data) 
//                      that can process different data types. One has to use the
//                      value returned by the RegisterFunction to identify the
//                      work type when providing new workload to the task manager.
//______________________________________________________________________________

class TThread;

class GeantTaskManager : public TNamed {
private:
   Int_t         fNactive;             // Number of active workers
   Bool_t        fToClean;             // Flag for requested cleanup
   vector<TThread*> fThreads;          // vector of threads
   vector<ExecuteTaskFunc_t>
                 exec_task;            // Vector of tasks of type void Function(void *data)
   tcqueue       fQueue;               // The work queue

   Bool_t        HasTask(ExecuteTaskFunc_t task) const;
   
public:
   GeantTaskManager(const char *name);
   virtual ~GeantTaskManager() {}
   
   Int_t         AddWorkers(Int_t nworkers);
   void          RemoveWorkers(Int_t nworkers, Bool_t sync=true);
   void          StopWorkers(Bool_t sync=true);
   
   Int_t         GetNtasks() const {return exec_task.size();}
   Int_t         GetNthreads() const {return fThreads.size()}
   Int_t         GetNactive() const {return fNworkers;}
   static void   ProcessLoop(void *arg);
   static void   MarkForRemoval(void *arg);

   // Main user methods
   // Register new task. Return task index.
   Int_t         RegisterTask(ExecuteTaskFunc_t task);
   // Add new data chunk to the processing queue. The task index to process the 
   // data must match the return value by RegisterTask
   void          ToProcess(void *data, Int_t task_id);
   

      
   ClassDef(GeantTaskManager, 0)      // A concurrent task manager
};
#endif
