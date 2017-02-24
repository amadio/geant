#include "DiscreteProcHandler.h"

#include "GeantTaskData.h"
#include "PhysicsProcessOld.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
DiscreteProcHandler::DiscreteProcHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
DiscreteProcHandler::~DiscreteProcHandler()
{
// Destructor
}  


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void DiscreteProcHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Invoke scalar post step processes

// ### this is to be refactored ###
      // Discrete processes only
/*
      nphys = output.SortByLimitingDiscreteProcess(); // only those that kPhysics and not continous limit
      if (nphys) {
        // reset number-of-interaction-legth-left
        for (auto itr = 0; itr < nphys; ++itr)
          output.fNintLenV[itr] = -1.0;

        // Do post step actions for particles suffering a given process.
        // Surviving particles are added to the output array
#ifdef USE_REAL_PHYSICS
        propagator->GetPhysicsInterface()->PostStepAction(mat, nphys, output, ntotnext, td);
#else
        // first: sample target and type of interaction for each primary tracks
        propagator->Process()->PostStepTypeOfIntrActSampling(mat, nphys, output, td);

//
// TODO: vectorized final state sampling can be inserted here through
//       a proper interface code that will:
//         1. select the necessary primary tracks based on the alrady
//            sampled interaction type
//         2. take all the member of the selected primary tracks that
//            necessary for sampling the final states
//         3. call the appropriate vector physics code that will
//            perform the physics interaction itself
//         4. update properties of the selected primary tracks,
//            insert secondary tracks into the track vector and set
//            'ntotnext' to the value of the number of secondary tracks
//            inserted to the track vector
//
  #ifdef USE_VECPHYS
        propagator->fVectorPhysicsProcess->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
  #endif
        // second: sample final states (based on the inf. regarding sampled
        //         target and type of interaction above), insert them into
        //         the track vector, update primary tracks;
        propagator->Process()->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
#endif
      }
*/


  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void DiscreteProcHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector post step handling.

  TrackVec_t &tracks = input.Tracks();
  // Copy tracks to output
#ifndef VECCORE_CUDA
  std::copy(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks) output->AddTrack(track);
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
