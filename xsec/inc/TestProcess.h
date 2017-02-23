#ifndef GEANT_TESTPROCESS
#define GEANT_TESTPROCESS

// Fast sim comments...

#include "Geant/Config.h"

#ifndef GEANT_PHYSICSPROCESS
#include "PhysicsProcessOld.h"
#endif

#include "base/Global.h"
#include "Geant/Typedefs.h"

#include "GeantFwd.h"


class TestProcess : public Geant::PhysicsProcessOld {
  public:

    using GeantTrack = Geant::GeantTrack;
    using GeantTrack_v = Geant::GeantTrack_v;
    using GeantTaskData = Geant::GeantTaskData;
    using TrackVec_t = Geant::TrackVec_t;

    TestProcess();
    virtual ~TestProcess() {}

    virtual void Initialize();
    virtual void ComputeIntLen( Material_t *mat, int ntracks, GeantTrack_v &tracks,
                                GeantTaskData *td );

    // dummy method: PostStep has been splitted up into two parts (see below)
    virtual void PostStep( Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks */, int & /*nout*/,
                           GeantTaskData * /*td*/ ) {}

    // sampling: target atom and type of the interaction for each primary tracks;
    //           all inf. regarding sampling output is stored in the tracks
    virtual void PostStepTypeOfIntrActSampling( Material_t *mat, int ntracks, GeantTrack_v &tracks,
                                                GeantTaskData *td );

    // sampling final states for each primary tracks based on target atom and
    // interaction type sampled by PostStepTypeOfIntrActSampling;
    // updating primary track properties and inserting secondary tracks;
    // number of inserted secondary tracks will be stored in nout at termination
    virtual void PostStepFinalStateSampling( Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                             GeantTaskData *td );

    virtual void AtRest( int /*ntracks*/, GeantTrack_v &/* tracks */, int &/* nout */, GeantTaskData */* td */ ) {};

    virtual void Eloss( Material_t */* mat */, int /* ntracks */, GeantTrack_v &/* tracks */, int &/* nout */, GeantTaskData */* td */ );

    virtual void ApplyMsc( Material_t */* mat */, int /* ntracks */, GeantTrack_v &/* tracks */, GeantTaskData */* td */ ) {};

//=== N E W   I N T E R F A C E S ===//
  VECCORE_ATT_HOST_DEVICE
  virtual void ComputeIntLen(Material_t *, TrackVec_t &, GeantTaskData *) {}
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepTypeOfIntrActSampling(Material_t *, TrackVec_t &, GeantTaskData *) {}
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepFinalStateSampling(Material_t *, TrackVec_t &, int &, GeantTaskData *) {}

  private:

    TestProcess( const TestProcess& );             // no imp.
    TestProcess& operator=( const TestProcess& );  // no imp.

};

#endif
