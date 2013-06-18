#ifndef GEANT_TASKBROKER
#define GEANT_TASKBROKER

class GeantTrack;

class TaskBroker
{
public:
   virtual void runTask(int threadid, int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin) = 0;
};


#endif // GEANT_TASKBROKER

