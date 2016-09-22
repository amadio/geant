// Implementation - for source file

//  This class interfaces between the GeantV Propagator and
//    the Compton model GUComptonKleinNishina
//
// Authors:  J. Apostolakis, M. Novak
#include "GVComptonProcess.h"

#include <iostream>
#include <math.h>
// Vector physics model related
#include "GUTrack.h"

// GeantV prototype related
#include "globals.h"
#include "GeantTrack.h"
#include "GeantTrackVec.h"
#include "GeantPropagator.h"

// Xsec library realted
#include "TPartIndex.h"
#include "base/PhysicalConstants.h" // from vecphys
#include "base/VecPhys.h" // for backend settings

//------------------------------------------------------------------------------
FQUALIFIER
GVComptonProcess::GVComptonProcess()
    : fProcessId(-1), fEnergyLimit(0.), fVComptonModel(0), fTargetElements(0), fParentTrackIndices(0),
      fPrimaryTracks(0), fSecondaryTracks(0) {}

//------------------------------------------------------------------------------
FQUALIFIER
GVComptonProcess::GVComptonProcess(int processId, double energyLimit, int numtracks)
    : fProcessId(processId), fEnergyLimit(energyLimit), fPrimaryTracks(0), fSecondaryTracks(0) {
  fPrimaryTracks = new GUTrack_v();
  fSecondaryTracks = new GUTrack_v();
  Allocator(numtracks);                         // allocate numtracks GUTrack_v for primaries
                                                // and secondaries and so on; will be increased if needed
  fVComptonModel = new vecphys::ComptonKleinNishina(0,-1); // create the vector physics model

  std::cout << "------------ process index: = " << fProcessId << std::endl;
}

//------------------------------------------------------------------------------
FQUALIFIER
GVComptonProcess::~GVComptonProcess() {
  Deallocator();
  delete fPrimaryTracks;
  delete fSecondaryTracks;
  delete fVComptonModel;
}

//------------------------------------------------------------------------------
FQUALIFIER
int GVComptonProcess::ApplyPostStepProcess(GeantTrack_v &gTrackV, int numtracks, int tid) {
  FilterPrimaryTracks(gTrackV, numtracks);
  if (fPrimaryTracks->numTracks == 0) // there is no track with Compton -> return
    return 0;
  PerformInteraction();
  return WriteBackTracks(gTrackV, tid);
}

//------------------------------------------------------------------------------
void GVComptonProcess::Allocator(int size) {
  if (size < 1)
    size = 1;
  fTargetElements = new int[size];
  fParentTrackIndices = new int[size];
  GUTrackAllocator(*fPrimaryTracks, size);
  GUTrackAllocator(*fSecondaryTracks, size);
}

//------------------------------------------------------------------------------
void GVComptonProcess::Deallocator() {
  if (fTargetElements)
    delete fTargetElements;
  fTargetElements = 0;
  if (fParentTrackIndices)
    delete fParentTrackIndices;
  fParentTrackIndices = 0;
  if (fPrimaryTracks)
    GUTrackDeallocator(*fPrimaryTracks);
  if (fSecondaryTracks)
    GUTrackDeallocator(*fSecondaryTracks);
}

//------------------------------------------------------------------------------
void GVComptonProcess::GUTrackAllocator(GUTrack_v &gutrack_v, int size) {
  if (gutrack_v.capacity > 0)
    GUTrackDeallocator(gutrack_v);

  // keep only those members that are really necessary
  gutrack_v.capacity = size; // maximum number of tracks that can be stored
  gutrack_v.numTracks = 0;   // real number of tracks stored
  //   gutrack_v.status        = new int[size];  // status of the tarcks stored: (we don't need this now) ???
  //   gutrack_v.particleType  = new int[size];  // type of the particles stored: what is this exactly ?
  //   gutrack_v.id            = new int[size];  // ??
  gutrack_v.parentId = new int[size]; // corresponding parent index in another GUTrack_v
  //   gutrack_v.proc          = new int[size];  // process index (we don't need this now) ???
  //   gutrack_v.x             = new double[size];   // (x,y,z) position (we don't need this now) ???
  //   gutrack_v.y             = new double[size];
  //   gutrack_v.z             = new double[size];
  gutrack_v.px = new double[size]; // momentum (px,py,pz)
  gutrack_v.py = new double[size];
  gutrack_v.pz = new double[size];
  gutrack_v.E = new double[size]; // total energy (kinetic E would be good as well)
  //   gutrack_v.q             = new double[size];   // charge (do we need this now ???)
  //   gutrack_v.s             = new double[size];   // current step length (do we need this now ???)
}

//------------------------------------------------------------------------------
void GVComptonProcess::GUTrackDeallocator(GUTrack_v &gutrack_v) {
  // keep only those members that are really necessary
  //   delete gutrack_v.status;
  //   delete gutrack_v.particleType;
  //   delete gutrack_v.id;
  delete gutrack_v.parentId;
  //   delete gutrack_v.proc;
  //   delete gutrack_v.x;
  //   delete gutrack_v.y;
  //   delete gutrack_v.z;
  delete gutrack_v.px;
  delete gutrack_v.py;
  delete gutrack_v.pz;
  delete gutrack_v.E;
  //   delete gutrack_v.q;
  //   delete gutrack_v.s;

  gutrack_v.capacity = 0;
}

//------------------------------------------------------------------------------
FQUALIFIER
void GVComptonProcess::FilterPrimaryTracks(GeantTrack_v &gTrackV, int numtracks) {
  // count number of tracks in GeantTrack_v with selected process = Compton
  int numRealTracks = 0;
  for (int i = 0; i < numtracks; ++i)
    if (gTrackV.fProcessV[i] == fProcessId)
      ++numRealTracks;

  // if no compton in this GeantTrack_v return immediately after
  // setting fPrimaryTracks->numTracks to 0
  if (numRealTracks < 1) {
    fPrimaryTracks->numTracks = 0;
    return;
  }

  // check if we have enough space: if not make sure that we have
  if (numRealTracks > fPrimaryTracks->capacity) {
    Deallocator();
    Allocator(numRealTracks);
  }

  // form the input GUTrack_v with the corresponding primary track members
  fPrimaryTracks->numTracks = numRealTracks;
  int j = 0;
  for (int i = 0; i < numtracks && j < numRealTracks; ++i) {
    if (gTrackV.fProcessV[i] == fProcessId) {               // this track suffered Compton: get it
      fParentTrackIndices[j] = i;                           // store index of this track in GeantTrack_v
                                                            // because we will need it to update
      double momentum = gTrackV.fPV[i];                     // total momentum
      fPrimaryTracks->px[j] = momentum * gTrackV.fXdirV[i]; // 3-momentum
      fPrimaryTracks->py[j] = momentum * gTrackV.fYdirV[i];
      fPrimaryTracks->pz[j] = momentum * gTrackV.fZdirV[i];
      fPrimaryTracks->E[j] = gTrackV.fEV[i]; // total energy

      fTargetElements[j] = gTrackV.fEindexV[i]; // Z of the target atom
      ++j;
    }
  }
}

//------------------------------------------------------------------------------
FQUALIFIER
void GVComptonProcess::PerformInteraction() {
  // call the vector physics model and perform the physics intarction:
  //   1. fPrimaryTracks is a GUTrack_v* that contains fPrimaryTracks->numTracks
  //      number of primary tracks
  //   2. fTargetElements is an int array with the corresponding Z of the target
  //      atoms
  //   3. fSecondaryTracks is a GUTrack_v* that can store fSecondaryTracks->capacity
  //      number of GUTrack-s (capacity >= fPrimaryTracks->numTracks now); the
  //      final number of secondary tracks produced by the vector physics model
  //      must be set in fSecondaryTracks->numTracks bythe model
  fVComptonModel->Interact<vecphys::VectorBackend>(*fPrimaryTracks, fTargetElements, *fSecondaryTracks);
}

//------------------------------------------------------------------------------
FQUALIFIER
int GVComptonProcess::WriteBackTracks(GeantTrack_v &gTrackV, int tid) {
  // 1. update primary tracks in GeantTrack_v to their post-interaction state
  int numPrimaryTracks = fPrimaryTracks->numTracks; // number of primary tracks has been used
  for (int ip = 0; ip < numPrimaryTracks; ++ip) {
    int indxP = fParentTrackIndices[ip]; // primary track indices in GeantTrack_v

    // set the selected process value for each Compton tracks in GeantTrack_v
    // to -1 i.e. no intercation: tabulated final states won't be sampled for
    // these tracks
    gTrackV.fProcessV[indxP] = -1;

    // check the post-interaction kinetic energy of the primary
    double kinE = fPrimaryTracks->E[ip] - gTrackV.fMassV[indxP]; // should be [GeV]
    if (kinE > fEnergyLimit) {                    // above tracking cut -> survived the phyisics intercation
      gTrackV.fEV[indxP] = fPrimaryTracks->E[ip]; // update total E [GeV]
      double totalP = std::sqrt(kinE * (kinE + 2.0 * gTrackV.fMassV[indxP])); // total P in [GeV]
      double invTotalP = 1.0 / totalP;
      gTrackV.fPV[indxP] = totalP; // update total P
      // assume that the physics model has already rotated the track (in a vectorized way);
      gTrackV.fXdirV[indxP] = fPrimaryTracks->px[ip] * invTotalP; // update x-dir
      gTrackV.fYdirV[indxP] = fPrimaryTracks->py[ip] * invTotalP; // update y-dir
      gTrackV.fZdirV[indxP] = fPrimaryTracks->pz[ip] * invTotalP; // update z-dir
    } else {                                                      // apply tracking cut:
      gTrackV.fEdepV[indxP] += kinE;                              // put ekin to energy depo
      gTrackV.fStatusV[indxP] = kKilled;                          // kill the primary track
      // should make sure that at-rest process is invoked if needed (not no ofc.)
    }
  }

  // 2. go for the secondaries
  int numSecondaries = fSecondaryTracks->numTracks;
  if (numSecondaries < 1) // if there are no secondaries return immediately
    return 0;

  // get the GeantPropagator
  GeantPropagator *propagator = GeantPropagator::Instance();

  // A process with the same secondary can define its base properties here
  // like Compton: the only possible secondary is e-
  const int secPDG = 11;                                    // e- PDG code
  //  const int secGVcode = TPartIndex::I()->PartIndex(secPDG); // e- GV code
  // the rest mass and charge of the secondary particle in a general way
  const double secMass = vecphys::electron_mass_c2;
  const double secCharge = vecphys::electron_charge;

  int numInsertedTracks = 0;
  for (int isec = 0; isec < numSecondaries; ++isec) {
    // primary track indices in GeantTrack_v
    int indxP = fParentTrackIndices[fSecondaryTracks->parentId[isec]];
    // check the kinetic energy of the secondary
    double kinE = fSecondaryTracks->E[isec] - secMass; // should be [GeV]
    if (kinE > fEnergyLimit) {                         // above tracking cut -> insert into GeantTrack_v
      // get a temporary GeantTrack from the propagator
      GeantTrack &gTrack = propagator->GetTempTrack(tid);

      // set some general properties: initial values or same as parent
      SetGeantTrack(gTrack, gTrackV, indxP);

      // set additianl members of gTrack
      double secPtot = std::sqrt(kinE * (kinE + 2. * secMass)); // total P [GeV]
      //      double invSecPtot = 1. / secPtot;                         // inv. total P [1/GeV]
      gTrack.fPDG = secPDG;                                     // PDG code
      gTrack.fGVcode = TPartIndex::I()->PartIndex(secPDG);      // corresponding GV code
      gTrack.fCharge = secCharge;                               // charge
      gTrack.fMass = secMass;                                   // rest mass [GeV]
      gTrack.fXdir = fSecondaryTracks->px[isec] * secPtot;      // direction (x,y,z)
      gTrack.fYdir = fSecondaryTracks->py[isec] * secPtot;
      gTrack.fZdir = fSecondaryTracks->pz[isec] * secPtot;
      gTrack.fP = secPtot;                   // total momentum
      gTrack.fE = fSecondaryTracks->E[isec]; // total energy

      // insert the new track
      propagator->AddTrack(gTrack);
      gTrackV.AddTrack(gTrack);

      // increase number of inserted tracks
      ++numInsertedTracks;
    } else { // apply tracking cut:
      // do not insert this track into GeantTrack_v
      gTrackV.fEdepV[indxP] += kinE; // put ekin to energy depo
      // should make sure that at-rest process is invoked if needed (not no ofc.)
    }
  }
  return numInsertedTracks;
}

//------------------------------------------------------------------------------
void GVComptonProcess::SetGeantTrack(GeantTrack &left, GeantTrack_v &right, int ip) {
  left.fEvent = right.fEventV[ip];   // same as parent
  left.fEvslot = right.fEvslotV[ip]; // same as parent
                                     //    left.fPDG      = secPDG;                              // will be set
                                     //    left.fGVcode   = TPartIndex::I()->PartIndex(secPDG);  // will be set
  left.fEindex = -1;                 // init
                                     //    left.fCharge   = secCharge;                           // will be set
  left.fProcess = -1;                // init
  //  left.fIzero = 0;                   // init
  left.fNsteps = 0;                  // init
  left.fStatus = kNew;               // new track
                                     //    left.fMass     = secMass;                             // will be set
  left.fXpos = right.fXposV[ip];     // position (same as parent)
  left.fYpos = right.fYposV[ip];
  left.fZpos = right.fZposV[ip];
  //    left.fXdir     = px*inv_Ptot;                         // dir. will be set
  //    left.fYdir     = py*inv_Ptot;
  //    left.fZdir     = pz*inv_Ptot;
  //    left.fP        = secPtot;                             // will be set
  //    left.fE        = secEtot;                             // will be set
  left.fEdep = 0.;                     // init
  left.fPstep = 0.;                    // init
  left.fStep = 0.;                     // init
  left.fSnext = 0.;                    // init
  left.fSafety = right.fSafetyV[ip];   // init to (same as parent)
  //  left.fFrombdr = right.fFrombdrV[ip]; // init to (same as parent)
  left.fPending = false;              // init
  *left.fPath = *right.fPathV[ip];     // init
  *left.fNextpath = *right.fPathV[ip]; // init
}
