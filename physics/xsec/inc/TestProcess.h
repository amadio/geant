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

    void Initialize() override;
    void ComputeIntLen( Material_t *mat, int ntracks, GeantTrack_v &tracks,
                        GeantTaskData *td ) override;

    // sampling: target atom and type of the interaction for each primary tracks;
    //           all inf. regarding sampling output is stored in the tracks
    void PostStepTypeOfIntrActSampling( Material_t *mat, int ntracks, GeantTrack_v &tracks,
                                        GeantTaskData *td ) override;

    // sampling final states for each primary tracks based on target atom and
    // interaction type sampled by PostStepTypeOfIntrActSampling;
    // updating primary track properties and inserting secondary tracks;
    // number of inserted secondary tracks will be stored in nout at termination
    void PostStepFinalStateSampling( Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                     GeantTaskData *td ) override;

    void AtRest( int /*ntracks*/, GeantTrack_v &/* tracks */, int &/* nout */, GeantTaskData */* td */ ) override {};

    void Eloss( Material_t */* mat */, int /* ntracks */, GeantTrack_v &/* tracks */, int &/* nout */, GeantTaskData */* td */ ) override;

    void ApplyMsc( Material_t */* mat */, int /* ntracks */, GeantTrack_v &/* tracks */, GeantTaskData */* td */ ) override {};

//=== N E W   I N T E R F A C E S ===//
  VECCORE_ATT_HOST_DEVICE
  virtual void ComputeIntLen(TrackVec_t &, GeantTaskData *) override {}
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepTypeOfIntrActSampling(TrackVec_t &, GeantTaskData *) override {}
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepFinalStateSampling(TrackVec_t &, int &, TrackVec_t &, GeantTaskData *) override {}

  private:

    TestProcess( const TestProcess& ) = delete;
    TestProcess& operator=( const TestProcess& ) = delete;

};

#endif
