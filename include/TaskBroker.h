#ifndef GEANT_TASKBROKER
#define GEANT_TASKBROKER

class GeantTrack;

class TaskBroker
{
protected:
   struct StreamHelper {};
public:
   typedef StreamHelper *Stream;
   virtual Stream GetNextStream() = 0;

   virtual void runTask(int threadid, int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin) = 0;

   virtual Stream launchTask(bool wait = false);
   virtual void waitForTasks();
   
   virtual long GetTotalWork();
};


#endif // GEANT_TASKBROKER

