#include "TTabPhysMgr.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#include "volumes/Particle.h"
using vecgeom::Particle;
#else
#include "TGeoManager.h"
#include "TGeoBranchArray.h"
#include "TGeoExtension.h"
#endif

#include "GeantTrack.h"

#include "globals.h"
#include "GeantPropagator.h"
#include "GeantTaskData.h"

#include "TRandom.h"
#include "TBits.h"
#include "TStopwatch.h"
#include "TError.h"
#include "TFile.h"
#include "TList.h"
#include "TSystem.h"

#include "TPartIndex.h"
#include "TEXsec.h"
#include "TMXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"

#include <iostream>

#include "base/RNG.h"
#include "Geant/Math.h"
using vecgeom::kPi;
using vecgeom::kTwoPi;

ClassImp(TTabPhysMgr)

    void loadvecgeomgeometry(GeantPropagator *);
TTabPhysMgr *TTabPhysMgr::fgInstance = 0;

//______________________________________________________________________________
TTabPhysMgr *TTabPhysMgr::Instance(const char *xsecfilename, const char *finalsfilename) {
  // Access to instance of TTabPhysMgr
  if (fgInstance)
    return fgInstance;
  if (!(xsecfilename && finalsfilename)) {
    ::Error("TTabPhysMgr::Instance", "Create TTabPhysMgr instance providing xsec files");
    return 0;
  }
  fgInstance = new TTabPhysMgr(xsecfilename, finalsfilename);
  return fgInstance;
}

//______________________________________________________________________________
TTabPhysMgr::~TTabPhysMgr() {
  // Destructor
  delete[] fMatXsec;
  delete[] fElemXsec;
  delete[] fElemFstate;
  delete fDecay;
  delete fHasNCaptureAtRest;
  fgInstance = 0;
}

//______________________________________________________________________________
TTabPhysMgr::TTabPhysMgr()
    : fNelements(0), fNmaterials(0), fElemXsec(0), fElemFstate(0), fMatXsec(0), fDecay(0),
#ifndef USE_VECGEOM_NAVIGATOR
      fGeom(0),
#endif
      fHasNCaptureAtRest(0) {
  // Dummy ctor.
  fgInstance = this;
}

//______________________________________________________________________________
TTabPhysMgr::TTabPhysMgr(const char *xsecfilename, const char *finalsfilename)
    : fNelements(0), fNmaterials(0), fElemXsec(0), fElemFstate(0), fMatXsec(0), fDecay(0),
#ifndef USE_VECGEOM_NAVIGATOR
      fGeom(0),
#endif
      fHasNCaptureAtRest(0) {

  fgInstance = this;
  TStopwatch timer;
  timer.Start();
// Load elements from geometry
#ifdef USE_VECGEOM_NAVIGATOR
  loadvecgeomgeometry(GeantPropagator::Instance());
  vector<vecgeom::Material *> matlist = vecgeom::Material::GetMaterials();
#else
  fGeom = gGeoManager;
  if (!fGeom)
    Fatal("TTabPhysMgr", "No geometry");
  TList *matlist = (TList *)fGeom->GetListOfMaterials();
  TIter next(matlist);
#endif
  Material_t *mat = 0;

  // Open xsec_FTFP_BERT.root file (or other phys.lists)
  TFile *fxsec = TFile::Open(xsecfilename);
  if (!fxsec) {
    Fatal("TTabPhysMgr", "Cannot open %s", xsecfilename);
  }
  fxsec->Get("PartIndex");
  // Open the fstate_FTFP_BERT.root file (or other phys.lists)
  TFile *fstate = TFile::Open(finalsfilename);
  if (!fstate) {
    Fatal("TTabPhysMgr", "Cannot open %s", finalsfilename);
  }

  // check version of the data files
  if (fgVersion != TPartIndex::I()->Version()) {
    std::cerr << "\n\n*************************************************************\n"
              << "  ---------------------------ERROR-----------------------------\n"
              << "    Your xsec_*.root and fstate_*.root data files at           \n"
              << "    -> " << xsecfilename << "\n"
              << "    -> " << finalsfilename << "\n"
              << "    Version is       : " << TPartIndex::I()->VersionMajor() << "." << TPartIndex::I()->VersionMinor()
              << "." << TPartIndex::I()->VersionSub() << "\n"
              << "    Required version : " << GetVersion() << "\n"
              << "    Update your xsec_*.root and fstate_*.root data files !     "
              << "\n*************************************************************\n\n";
    exit(EXIT_FAILURE);
  }

  // get the decay table from the final state file
  fDecay = (TPDecay *)fstate->Get("DecayTable");

#ifdef USE_VECGEOM_NAVIGATOR
  printf("#materials:= %lu \n", matlist.size());
#else
  // INFO: print number of materials in the current GeoManager
  printf("#materials:= %d \n", matlist->GetSize());
#endif

  // First loop on all materials to mark used elements
  TBits elements(NELEM);
#ifdef USE_VECGEOM_NAVIGATOR
  for (int i = 0; i < matlist.size(); ++i) {
    mat = matlist[i];
#else
  while ((mat = (Material_t *)next())) {
#endif
    std::cout << mat->GetName() << "used " << mat->IsUsed() << " Z " << mat->GetZ() << std::endl;
    if (!mat->IsUsed() || mat->GetZ() < 1.)
      continue;
    fNmaterials++;
    int nelem = mat->GetNelements();
    // Check if we are on the safe side; should exit otherwise
    if (nelem > MAXNELEMENTS) {
      Fatal("TTabPhysMgr", "Number of elements in %s is %d > TTabPhysMgr::MAXNELEMENTS=%d\n", mat->GetName(), nelem,
            MAXNELEMENTS);
    }
    for (int iel = 0; iel < nelem; ++iel) {
      double ad;
      double zd;
      double wd;
      mat->GetElementProp(ad, zd, wd, iel);
      if (zd < 1 || zd > NELEM) {
        Fatal("TTabPhysMgr", "In material %s found element with z=%d > NELEM=%d", mat->GetName(), (int)zd, NELEM);
      }
      elements.SetBitNumber(zd);
    }
  }
  fNelements = elements.CountBits();
  fElemXsec = new TEXsec *[NELEM];
  fElemFstate = new TEFstate *[NELEM];
  fMatXsec = new TMXsec *[fNmaterials];
  printf("Reading xsec data and final states for %d elements in %d materials\n", fNelements, fNmaterials);
  // Loop elements and load corresponding xsec and final states
  int zel = elements.FirstSetBit();
  int nbits = elements.GetNbits();
  TEXsec *exsec;
  TEFstate *estate;
  // Load elements xsec data in memory
  ProcInfo_t procInfo1, procInfo2;
  gSystem->GetProcInfo(&procInfo1);
  while (zel < nbits) {
    exsec = TEXsec::GetElement(zel, 0, fxsec);
    fElemXsec[zel] = exsec;
    fElemXsec[zel]->SetIndex(zel); // for quick access to the corresponding fstate
    estate = TEFstate::GetElement(zel, 0, fstate);
    fElemFstate[zel] = estate;
    printf("   loaded xsec data and states for: %s\n", TPartIndex::I()->EleSymb(zel));
    zel = elements.FirstSetBit(zel + 1);
    // init : does the particle have nuclear cpature at rest? array
    if (!fHasNCaptureAtRest) {
      int numParticles = TPartIndex::I()->NPart();
      fHasNCaptureAtRest = new bool[numParticles];
      for (int ip = 0; ip < numParticles; ++ip)
        fHasNCaptureAtRest[ip] = estate->HasRestCapture(ip);
    }
  }
  gSystem->GetProcInfo(&procInfo2);
  long mem = (procInfo2.fMemResident - procInfo1.fMemResident) / 1024;
  fxsec->Close();
  fstate->Close();
  // xsec and states now in memory
  // Go through all materials in the geometry and form the associated TMXsec
  // objects.
  int *z = new int[MAXNELEMENTS];
  int *a = new int[MAXNELEMENTS];
  float *w = new float[MAXNELEMENTS];
  fNmaterials = 0;
#ifdef USE_VECGEOM_NAVIGATOR
  for (int i = 0; i < matlist.size(); ++i) {
    mat = matlist[i];
#else
  next.Reset();
  while ((mat = (Material_t *)next())) {
#endif
    //    std::cout << __FILE__ << "::" << __func__ << "::Loading xsec for " << mat->GetName() << std::endl;
    if (!mat->IsUsed())
      continue;
    int nelem = mat->GetNelements();
    // loop over the elements of the current material in order to obtain the
    // z, a, w, arrays of the elements of this material
    double ad;
    double zd;
    double wd;
    for (int iel = 0; iel < nelem; ++iel) {
      mat->GetElementProp(ad, zd, wd, iel);
      a[iel] = ad;
      z[iel] = zd;
      w[iel] = wd;
    }
    if (nelem == 0) {
      mat->Dump();
      Fatal("TTabPhysMgr", "The material (%s) seems to have no elements", mat->GetName());
    }
    // Construct the TMXsec object that corresponds to the current material
    TMXsec *mxs = new TMXsec(mat->GetName(), mat->GetName(), z, a, w, nelem, mat->GetDensity(), kTRUE, fDecay);
    fMatXsec[fNmaterials++] = mxs;
// Connect to Material
#ifdef USE_VECGEOM_NAVIGATOR
    mat->SetXsecPtr(static_cast<void *>(mxs));
#else
    mat->SetFWExtension(new TGeoRCExtension(new TOMXsec(mxs)));
#endif
  } // End of while
  delete[] z;
  delete[] a;
  delete[] w;

  // After setting up all the necessary TMXsec objects we have the arra of the
  // loaded elemental TEXsec object pointers in: static TEXsec::TEXsec *fElements[NELEM]
  // Since the static TEXsec *fElements[NELEM] is private in TEXsec, I added a getter:
  // IN TEXsec.h:
  // static TEXsec** GetElements(){ return fElements; }
  // that I will use here to set TTabPhysMgr::fElemXsec )
  // fElemXsec = TEXsec::GetElements();
  int nelements = TEXsec::NLdElems();
  if (nelements != fNelements)
    Error("TTabPhysMgr", "Number of elements not matching");

  // INFO: print some info for checking
  printf("number of materials in fMatXsec[]:= %d\n", fNmaterials);
  for (int i = 0; i < fNmaterials; ++i)
    printf("   fMatXsec[%d]: %s\n", i, fMatXsec[i]->GetName());
  timer.Stop();
  printf("Memory taken by xsec and states: %ld [MB] loaded in: %g [sec]\n", mem, timer.CpuTime());
}

//______________________________________________________________________________
void TTabPhysMgr::TransformLF(int /*indref*/, GeantTrack_v & /*tracks*/, int /*nproducts*/, int /*indprod*/,
                              GeantTrack_v & /*output*/) {
  // Transform tracks taken from the final state from the local frame to the lab
  // frame (LF). Not clear what parameters to add yet.
  // Input: reference track (mother) described as vector container + index of ref track
  // Input: number of tracks in the final state, start index and vector container
  // Output: roto-boosted tracks in the output vector
}

// NOT ACTIVE NOW
//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
void TTabPhysMgr::ApplyMsc(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Compute MSC angle at the beginning of the step and apply it to the vector
  // of tracks.
  // Input: material index, number of tracks in the tracks vector to be used
  // Output: fXdirV, fYdirV, fZdirV modified in the track container for ntracks

  TMXsec *mxs = 0;
#ifndef GEANT_CUDA_DEVICE_BUILD
  if (mat)
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)mat->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)mat->GetFWExtension())->GetUserObject())->MXsec();
#endif
#else
  // NOTE: we need to get it from somewhere ....
  assert(mxs != 0);
#endif
  //   static int icnt=0;
  double msTheta;
  double msPhi;

#ifndef GEANT_CUDA_DEVICE_BUILD
  double *rndArray = td->fDblArray;
  td->fRndm->RndmArray(ntracks, rndArray);
#else
  double *rndArray = 0; // NOTE: we need to get it from somewhere ....
  VECGEOM_NAMESPACE::RNG::Instance().uniform_array(ntracks, rndArray, 0., 1.);
#endif

  //   double dir[3] = {0.,0.,0.};
  if (mxs) {
    for (int i = 0; i < ntracks; ++i) {
      msTheta = mxs->MS(tracks.fGVcodeV[i], tracks.fEV[i] - tracks.fMassV[i]);
      msPhi = 2. * kPi * rndArray[i];
      RotateTrack(tracks, i, msTheta, msPhi);
    }
    return;
  }
  // Mixed tracks in different volumes
  for (int i = 0; i < ntracks; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)tracks.GetMaterial(i)->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)tracks.GetMaterial(i)->GetFWExtension())->GetUserObject())->MXsec();
#endif
    msTheta = mxs->MS(tracks.fGVcodeV[i], tracks.fEV[i] - tracks.fMassV[i]);
    msPhi = 2. * kPi * rndArray[i];
    /*
          if (icnt<100 && mat->GetZ()>10) {
             Printf("theta=%g  phi=%g", msTheta*vecgeom::Materialh::RadToDeg(), msPhi*vecgeom::Materialh::RadToDeg());
             dir[0] = tracks.fXdirV[i];
             dir[1] = tracks.fYdirV[i];
             dir[2] = tracks.fZdirV[i];
          }
    */
    RotateTrack(tracks, i, msTheta, msPhi);
    /*
          if (icnt<100 && mat->GetZ()>10) {
             icnt++;
             double dot = dir[0]*tracks.fXdirV[i] + dir[1]*tracks.fYdirV[i] +dir[2]*tracks.fZdirV[i];
             double angle = vecgeom::Materialh::ACos(dot)*vecgeom::Materialh::RadToDeg();
             Printf("new angle=%g   delta=%g", angle,
       vecgeom::Materialh::Abs(angle-msTheta*vecgeom::Materialh::RadToDeg()));
          }
    */
  }
}

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
int TTabPhysMgr::Eloss(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Apply energy loss for the input material for ntracks in the vector of
  // tracks. Output: modified tracks.fEV array

  TMXsec *mxs = 0;
#ifndef GEANT_CUDA_DEVICE_BUILD
  if (mat)
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)mat->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)mat->GetFWExtension())->GetUserObject())->MXsec();
#endif
#else
  // NOTE: we need to get it from somewhere ....
  assert(mxs != 0);
#endif
  int nTotSecPart = 0; // total number of new tracks
  double energyLimit = gPropagator->fEmin;
  if (mxs) {
    mxs->Eloss(ntracks, tracks);
    // call atRest sampling for tracks that have been stopped by Eloss and has at-rest
    for (int i = 0; i < ntracks; ++i)
      if (tracks.fProcessV[i] == -2 && HasRestProcess(tracks.fGVcodeV[i]))
        GetRestFinStates(tracks.fGVcodeV[i], mxs, energyLimit, tracks, i, nTotSecPart, td);
    return nTotSecPart;
  }
  // Mixed tracks in different volumes
  for (int i = 0; i < ntracks; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)tracks.GetMaterial(i)->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)tracks.GetMaterial(i)->GetFWExtension())->GetUserObject())->MXsec();
#endif
    mxs->ElossSingle(i, tracks);
    // call atRest sampling for tracks that have been stopped by Eloss and has at-rest
    if (tracks.fProcessV[i] == -2 && HasRestProcess(tracks.fGVcodeV[i]))
      GetRestFinStates(tracks.fGVcodeV[i], mxs, energyLimit, tracks, i, nTotSecPart, td);
  }

  return nTotSecPart;
}

//______________________________________________________________________________
void TTabPhysMgr::ProposeStep(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Sample free flight/proposed step for the firts ntracks tracks and store them
  // in tracks.fPstepV

  TMXsec *mxs = 0;
#ifndef GEANT_CUDA_DEVICE_BUILD
  if (mat) {
#ifdef USE_VECGEOM_NAVIGATOR
    /*
    static std::mutex m;
    m.lock();
    std::cout << __FILE__ << "::" << __func__ <<"::mat:" << mat->GetName()
         << " vol:" << vecgeom::GeoManager::Instance().FindLogicalVolume(tracks.fVindexV[0])->GetName()
         << " xsecptr:" << mat->GetXsecPtr() << std::endl;
    m.unlock();
    */
    mxs = (TMXsec *)mat->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)mat->GetFWExtension())->GetUserObject())->MXsec();
#endif
  }
#endif
  if (mxs) {
    mxs->ProposeStep(ntracks, tracks, td);
    return;
  }
  // Mixed tracks in different volumes
  for (int i = 0; i < ntracks; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)tracks.GetMaterial(i)->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)tracks.GetMaterial(i)->GetFWExtension())->GetUserObject())->MXsec();
#endif
    mxs->ProposeStepSingle(i, tracks, td);
  }
}

// Implemented in a different way
//______________________________________________________________________________
int TTabPhysMgr::SampleDecay(int /*ntracks*/, GeantTrack_v & /*tracksin*/, GeantTrack_v & /*tracksout*/) {
  // Sample decay for the tracks in the input vector and push the resulting tracks in
  // the output vector. Change status of decayed tracks. Returns number of new tracks.
  return 0;
}

// smapling: target atom and type of the interaction for each primary tracks
//______________________________________________________________________________
void TTabPhysMgr::SampleTypeOfInteractions(int imat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  Material_t *mat = 0;
  TMXsec *mxs = 0;
  if (imat >= 0) {
#ifdef USE_VECGEOM_NAVIGATOR
    mat = vecgeom::Material::GetMaterials()[imat];
    mxs = (TMXsec *)mat->GetXsecPtr();
#else
    mat = (TGeoMaterial *)fGeom->GetListOfMaterials()->At(imat);
    mxs = ((TOMXsec *)((TGeoRCExtension *)mat->GetFWExtension())->GetUserObject())->MXsec();
#endif
  }

  // 1. sampling: a. decay or something else
  //             b. if else then what on what target?
  //   output of sampling is stored in the tracks
  if (mxs) {
    mxs->SampleInt(ntracks, tracks, td);
    return;
  }
  // Mixed tracks in different volumes
  for (int i = 0; i < ntracks; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)tracks.GetMaterial(i)->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)tracks.GetMaterial(i)->GetFWExtension())->GetUserObject())->MXsec();
#endif
    mxs->SampleSingleInt(i, tracks, td);
  }
}

//______________________________________________________________________________
int TTabPhysMgr::SampleFinalStates(int imat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  GeantPropagator *propagator = GeantPropagator::Instance();
  double energyLimit = propagator->fEmin;

  Material_t *mat = 0;
  TMXsec *mxs = 0;
  if (imat >= 0) {
#ifdef USE_VECGEOM_NAVIGATOR
    mat = vecgeom::Material::GetMaterials()[imat];
    mxs = (TMXsec *)mat->GetXsecPtr();
#else
    mat = (Material_t *)fGeom->GetListOfMaterials()->At(imat);
    mxs = ((TOMXsec *)((TGeoRCExtension *)mat->GetFWExtension())->GetUserObject())->MXsec();
#endif
  }

  // tid-based rng
  double *rndArray = td->fDblArray;
  td->fRndm->RndmArray(2 * ntracks, rndArray);

  int nTotSecPart = 0; // total number of secondary particles in tracks

  for (int t = 0; t < ntracks; ++t) {
    // if no interaction was selected for this track (because it doesn't have any)
    if (tracks.fProcessV[t] < 0)
      continue;

    int nSecPart = 0;     // number of secondary particles per reaction
    const int *pid = 0;   // GeantV particle codes [nSecPart]
    const float *mom = 0; // momentum vectors the secondaries [3*nSecPart]
    float ener = 0;       // energy at the fstate (Ekin of primary after the interc.)
    float kerma = 0;      // released energy
    float weight = 0;     // weight of the fstate (just a dummy parameter now)
    char isSurv = 0;      // is the primary survived the interaction
    int ebinindx = -1;    // energy bin index of the selected final state

    // Deal with mixed tracks case
    if (!mxs)
#ifdef USE_VECGEOM_NAVIGATOR
      mxs = (TMXsec *)tracks.GetMaterial(t)->GetXsecPtr();
#else
      mxs = ((TOMXsec *)((TGeoRCExtension *)tracks.GetMaterial(t)->GetFWExtension())->GetUserObject())->MXsec();
#endif

    // firts check the results of interaction sampling:
    if (tracks.fProcessV[t] == 3) {
      // decay : in-flight decay was selected
      // kill the primary tarck
      tracks.fStatusV[t] = kKilled;
      // sample in-flight decay final state
      SampleDecayInFlight(tracks.fGVcodeV[t], mxs, energyLimit, tracks, t, nTotSecPart, td);
      tracks.fPV[t] = 0.;
      tracks.fEV[t] = tracks.fMassV[t];
      continue;
    }

    // not decay but something else was selected
    double curPrimEkin = tracks.fEV[t] - tracks.fMassV[t];
    isSurv = fElemFstate[tracks.fEindexV[t]]->SampleReac(tracks.fGVcodeV[t], tracks.fProcessV[t], curPrimEkin, nSecPart,
                                                         weight, kerma, ener, pid, mom, ebinindx, rndArray[2 * t],
                                                         rndArray[2 * t + 1]);

    // it is the case of: pre-step energy sigma is not zero of this interaction
    //                    but post step is zero-> we don't have final state for
    //                    this interaction at the postStep energy bin-> do nothing
    // This can happen if interaction type is selected based on the pre-step energy
    // and there is some energy loss along the step (or something like this)
    if (isSurv && ener < 0) // let it go further as it is
      continue;

    // we should correct the kerma as well but we don't have enough information
    tracks.fEdepV[t] += kerma;

    // if we have secondaries from the current interaction
    if (nSecPart) {
      double oldXdir = tracks.fXdirV[t]; // old X direction of the primary
      double oldYdir = tracks.fYdirV[t]; // old Y direction of the primary
      double oldZdir = tracks.fZdirV[t]; // old Z direction of the primary
      int j = 0;

      // setting the final state correction factor (we scale only the 3-momentums)
      //-get mass of the primary
      double primMass = tracks.fMassV[t]; // mass [GeV]
      //-compute corFactor = P_current/P_original = Pz_current/Pz_original
      // (normaly a check would be good but not necessary: if(ebinindx<0 -> ...)
      double orgPrimEkin = (TPartIndex::I()->EGrid())[ebinindx];
      double corFactor =
          sqrt(curPrimEkin * (curPrimEkin + 2.0 * primMass) / (orgPrimEkin * (orgPrimEkin + 2.0 * primMass)));
      //-if corFactor is set here to 1.0 --> no correction of the final states
      // corFactor = 1.0;

      // check if we need to correct the post-interaction Ekin of the primary:
      // if the primary is survived and has non-zero Ekin --> compute its corrected Ekin
      double postEkinOfParimary = ener;
      if (isSurv && (postEkinOfParimary > 0.0)) { // survived and not stopped
        // get corrected 3-momentum of the post-interaction primary
        double px = mom[0];
        double py = mom[1];
        double pz = mom[2];
        px *= corFactor;
        py *= corFactor;
        pz *= corFactor;
        // compute corrected P^2 in [GeV^2]
        double postPrimP2 = px * px + py * py + pz * pz;
        // recompute post-interaction Ekin of the primary with corrected 3-momentum
        postEkinOfParimary = sqrt(postPrimP2 + primMass * primMass) - primMass;
      }

      if (postEkinOfParimary > energyLimit) { // survived even after the correction and the E-limit.
        // keep alive
        //        tracks.fStatusV[t] = kAlive;
        double px = mom[0];
        double py = mom[1];
        double pz = mom[2];
        px *= corFactor;
        py *= corFactor;
        pz *= corFactor;
        // compute corrected P^2 in [GeV^2]
        double postPrimP2 = px * px + py * py + pz * pz;
        // recompute post-interaction Ekin of the primary with corrected 3-momentum
        postEkinOfParimary = sqrt(postPrimP2 + primMass * primMass) - primMass;

        // update primary in tracks
        double secPtot = sqrt(postPrimP2);                // total P [GeV]
        double secEtot = postEkinOfParimary + tracks.fMassV[t]; // total energy in [GeV]
        tracks.fPV[t] = secPtot;                                  // momentum of this particle
        tracks.fEV[t] = secEtot;                                  // total E of this particle
        tracks.fXdirV[t] = px / secPtot;                          // dirx of this particle (before transform.)
        tracks.fYdirV[t] = py / secPtot;                          // diry of this particle (before transform.)
        tracks.fZdirV[t] = pz / secPtot;                          // dirz of this particle (before transform.)

        // Rotate parent track in tracks to original parent track's frame
        RotateNewTrack(oldXdir, oldYdir, oldZdir, tracks, t);
        // primary track is updated
      } else {
        // Primary particle energy is below tracking limit
        //-set status of primary in tracks to kKilled;
        tracks.fStatusV[t] = kKilled;
        tracks.fEdepV[t] += postEkinOfParimary;
        tracks.fPV[t] = 0.;
        tracks.fEV[t] = tracks.fMassV[t];
        // if the primary is stopped i.e. Ekin <= 0 then call at-rest if it has
        if (isSurv && postEkinOfParimary <= 0.0 && HasRestProcess(tracks.fGVcodeV[t]))
          GetRestFinStates(tracks.fGVcodeV[t], mxs, energyLimit, tracks, t, nTotSecPart, td);
      }

      if (isSurv)
        j = 1;
      // loop over the secondaries and put them into tracks if they good to track:
      // j=0 -> including stopped primary as well if isSurv = kTRUE;
      // j=1 -> skipp the primary in the list of secondaries (was already updated in tracks above)
      for (int i = j; i < nSecPart; ++i) {
        if (pid[i] >= TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
          int idummy = pid[i] - 1000000000;
          int Z = idummy / 10000.;
          int A = (idummy - Z * 10000) / 10.;
          double secMass = TPartIndex::I()->GetAprxNuclearMass(Z, A);
          double px = mom[3 * i];
          double py = mom[3 * i + 1];
          double pz = mom[3 * i + 2];
          px *= corFactor;
          py *= corFactor;
          pz *= corFactor;
          double secPtot2 = px * px + py * py + pz * pz; // total P^2 [GeV^2]
          tracks.fEdepV[t] += sqrt(secPtot2 + secMass * secMass) - secMass;
          continue;
        }
        int secPDG = TPartIndex::I()->PDG(pid[i]); // Geant V particle code -> particle PGD code
#ifdef USE_VECGEOM_NAVIGATOR
        const Particle *const &secPartPDG = &Particle::GetParticle(secPDG);
#else
        TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
#endif
        double secMass = secPartPDG->Mass();
        /*	static std::mutex m;
        m.lock();
        std::cout << __func__ << "::secMass: " << secMass << " secPDG: " << secPDG << " SecPartPDG:" << *secPartPDG <<
        std::endl;
        m.unlock();
        */
        double px = mom[3 * i];
        double py = mom[3 * i + 1];
        double pz = mom[3 * i + 2];
        px *= corFactor;
        py *= corFactor;
        pz *= corFactor;
        double secPtot2 = px * px + py * py + pz * pz;             // total P^2 [GeV^2]
        double secPtot = sqrt(secPtot2);                     // total P [GeV]
        double secEtot = sqrt(secPtot2 + secMass * secMass); // total energy in [GeV]
        double secEkin = secEtot - secMass;                        // kinetic energy in [GeV]
        // Ekin of the i-th secondary is higher than the threshold
        if (secEkin >= energyLimit) { // insert secondary into OUT tracks_v and rotate
          GeantTrack &gTrack = td->GetTrack();
          //          GeantTrack gTrack;
          // set the new track properties
          gTrack.fEvent = tracks.fEventV[t];
          gTrack.fEvslot = tracks.fEvslotV[t];
          //          gTrack.fParticle = nTotSecPart;          //index of this particle
          gTrack.fPDG = secPDG;    // PDG code of this particle
          gTrack.fGVcode = pid[i]; // GV index of this particle
          gTrack.fEindex = 0;
          gTrack.fCharge = secPartPDG->Charge() / 3.; // charge of this particle
          gTrack.fProcess = 0;
          gTrack.fVindex = tracks.fVindexV[t];
          gTrack.fNsteps = 0;
          //          gTrack.fSpecies  = 0;
          gTrack.fStatus = kNew;           // status of this particle
          gTrack.fMass = secMass;          // mass of this particle
          gTrack.fXpos = tracks.fXposV[t]; // rx of this particle (same as parent)
          gTrack.fYpos = tracks.fYposV[t]; // ry of this particle (same as parent)
          gTrack.fZpos = tracks.fZposV[t]; // rz of this particle (same as parent)
          gTrack.fXdir = px / secPtot;     // dirx of this particle (before transform.)
          gTrack.fYdir = py / secPtot;     // diry of this particle before transform.)
          gTrack.fZdir = pz / secPtot;     // dirz of this particle before transform.)
          gTrack.fP = secPtot;             // momentum of this particle
          gTrack.fE = secEtot;             // total E of this particle
          gTrack.fTime = tracks.fTimeV[t]; // global time
          gTrack.fEdep = 0.;
          gTrack.fPstep = 0.;
          gTrack.fStep = 0.;
          gTrack.fSnext = 0.;
          gTrack.fSafety = tracks.fSafetyV[t];
          gTrack.fFrombdr = tracks.fFrombdrV[t];
          gTrack.fPending = kFALSE;
          *gTrack.fPath = *tracks.fPathV[t];
          *gTrack.fNextpath = *tracks.fPathV[t];

          // Rotate new track to parent track's frame
          RotateNewTrack(oldXdir, oldYdir, oldZdir, gTrack);

          propagator->AddTrack(gTrack);
          tracks.AddTrack(gTrack);

          ++nTotSecPart;
        } else {                       // {secondary Ekin < energyLimit} -> kill this secondary
          tracks.fEdepV[t] += secEkin; // add the Ekin of this secondary to the energy depositon
          // is secEkin <=0 then call at-rest process if the sec. particle has any
          if (secEkin <= 0.0 && HasRestProcess(pid[i]))
            GetRestFinStates(pid[i], mxs, energyLimit, tracks, t, nTotSecPart, td);
        }
      }                             // end loop over the secondaries
    } else {                        // nSecPart = 0 i.e. there is no any secondaries -> primary was killed as well
      tracks.fStatusV[t] = kKilled; // set status of primary in tracks to kKilled;
      tracks.fPV[t] = 0;
      tracks.fEV[t] = tracks.fMassV[t];
    }
  } // end loop over tracks

  return nTotSecPart;
}

/*
//______________________________________________________________________________
int TTabPhysMgr::SampleInt(int imat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td)
{
// 0. ntracks contains particles with status of Alive
// 1.Sampling the element of the material for interaction based on the relative
// total X-secs of the elements; Sampling the type of the interaction (on the
// sampled element) based on the realtive total X-secs of the interactions ;
// OUT:-indices of the TEXsec* in fElemXsec, that correspond to the sampled
//      elements, will be in GeantTrack_v::fEindexV array; GeantTrack_v::fEindexV[i]
//      will be -1 if no reaction for i-th particle
//     -the GV reaction indices will be in GeantTrack_v::fProcessV array;
//      GeantTrack_v::fEindexV[i] will be -1 if no reaction for i-th particle
// 2.Sampling the finale states for the selected interaction and store the secondary
// tracks in tracks; only those traks go into tracks that can pass the energyLimit,
// Rest process final states will be sampled in case of those secondaries that
// stopped. So the status of tracks can be only kAlive in tracks
// 4.number of secondary tracks will be returned and original track status will
// be updated (if they have been killed)

   // # smapling: target atom and type of the interaction for each primary tracks
   SampleTypeOfInteractions(imat, ntracks, tracks, td);

   // # sampling final states for each primary tracks based on target atom and
   //    interaction type sampled in SampleTypeOfInteractionsInt;
   // # upadting primary track properties and inserting secondary tracks;
   // # return: number of inserted secondary tracks
   return SampleFinalStates(imat, ntracks, tracks, td);
}
*/

// Will be called only if the particle has decay or/and nuclear capture at-rest
//______________________________________________________________________________
// will be called recursively if necessary
GEANT_CUDA_DEVICE_CODE
void TTabPhysMgr::GetRestFinStates(int partindex, TMXsec *mxs, double energyLimit, GeantTrack_v &tracks,
                                   int iintrack, int &nTotSecPart, GeantTaskData *td) {
  // current track should have already been killed before calling
  const double mecc = 0.00051099906; // e- mass c2 in [GeV]
  double rndArray[3];
#ifndef GEANT_CUDA_DEVICE_BUILD
  td->fRndm->RndmArray(3, rndArray);
#else
  VECGEOM_NAMESPACE::RNG::Instance().uniform_array(3, rndArray, 0., 1.);
#endif

  int nSecPart = 0;     // number of secondary particles per reaction
  const int *pid = 0;   // GeantV particle codes [nSecPart]
  const float *mom = 0; // momentum vectors the secondaries [3*nSecPart]
  float ener = 0;       // energy at the fstate
  float kerma = 0;      // released energy
  float weight = 0;     // weight of the fstate (just a dummy parameter now)
  char isSurv = 0;      // is the primary survived the interaction

  // check if particle is e+ : e+ annihilation at rest if $E_{limit}< m_{e}c^{2}$
  if (partindex == TPartIndex::I()->GetSpecGVIndex(1)) {
    if (energyLimit < mecc) {
      double randDirZ = 1.0 - 2.0 * rndArray[0];
      double randSinTheta = sqrt(1.0 - randDirZ * randDirZ);
      double randPhi = 2.0 * rndArray[1] * kPi;
      double randDirX = randSinTheta * cos(randPhi);
      double randDirY = randSinTheta * sin(randPhi);

      // need to do it one-by-one
      // 1. gamma
      GeantTrack &gTrack1 = td->GetTrack();
      // set the new track properties: 2 gamma with m_{e}*c*c
      gTrack1.fEvent = tracks.fEventV[iintrack];
      gTrack1.fEvslot = tracks.fEvslotV[iintrack];
      //       gTrack.fParticle = nTotSecPart;          //index of this particle
      gTrack1.fPDG = 22;                                    // gamma PDG code
      gTrack1.fGVcode = TPartIndex::I()->GetSpecGVIndex(2); // gamma GV index
      gTrack1.fEindex = 0;
      gTrack1.fCharge = 0.; // charge
      gTrack1.fProcess = 0;
      gTrack1.fVindex = tracks.fVindexV[iintrack];
      gTrack1.fNsteps = 0;
      //       gTrack.fSpecies  = 0;
      gTrack1.fStatus = kNew;                  // status of this particle
      gTrack1.fMass = 0.;                      // mass of this particle
      gTrack1.fXpos = tracks.fXposV[iintrack]; // rx of this particle (same as parent)
      gTrack1.fYpos = tracks.fYposV[iintrack]; // ry of this particle (same as parent)
      gTrack1.fZpos = tracks.fZposV[iintrack]; // rz of this particle (same as parent)
      gTrack1.fXdir = randDirX;
      gTrack1.fYdir = randDirY;
      gTrack1.fZdir = randDirZ;
      gTrack1.fP = mecc;                       // momentum of this particle
      gTrack1.fE = mecc;                       // total E of this particle
      gTrack1.fTime = tracks.fTimeV[iintrack]; // total time of this particle
      gTrack1.fEdep = 0.;
      gTrack1.fPstep = 0.;
      gTrack1.fStep = 0.;
      gTrack1.fSnext = 0.;
      gTrack1.fSafety = tracks.fSafetyV[iintrack];
      gTrack1.fFrombdr = tracks.fFrombdrV[iintrack];
      gTrack1.fPending = kFALSE;
      *gTrack1.fPath = *tracks.fPathV[iintrack];
      *gTrack1.fNextpath = *tracks.fPathV[iintrack];

      gPropagator->AddTrack(gTrack1);
      tracks.AddTrack(gTrack1);

      // 2. gamma : everything is the same but the direction
      gTrack1.fXdir = -1. * randDirX;
      gTrack1.fYdir = -1. * randDirY;
      gTrack1.fZdir = -1. * randDirZ;

      gPropagator->AddTrack(gTrack1);
      tracks.AddTrack(gTrack1);

      nTotSecPart += 2;
      return;
    } else {
      return;
    }
  }
  // If the stopped particle doesn't have nuclear capture at-rest then decay it
  if (!fHasNCaptureAtRest[partindex]) {
    // Decay at-rest
    // sample final state for decay
    isSurv = fDecay->SampleDecay(partindex, nSecPart, pid, mom);
  } else {
    // It has nuclear capture at rest so invoke that
    // sample one element of the material
    TEFstate *elemfstate = fElemFstate[mxs->SampleElement(td)];
    // stample final state for nuclear capture at-rest
    isSurv = elemfstate->SampleRestCaptFstate(partindex, nSecPart, weight, kerma, ener, pid, mom, rndArray[0]);
  }

  double randDirX = 0;
  double randDirY = 0;
  double randDirZ = 1;
  double randSinTheta;
  double randPhi;

  // note: parent was already stopped because an at Rest process happend;
  //      -> primary is not in the list of secondaries

  // isSurv should always be kFALSE here because primary was stopped -> just a check
  if (isSurv)
    printf("A stopped particle survived its rest process in TTabPhysMgr::GetRestFinSates!\n");

  // for a random rotation
  if (nSecPart) {
    randDirZ = 1.0 - 2.0 * rndArray[1];
    randSinTheta = sqrt((1.0 - randDirZ) * (1.0 + randDirZ));
    randPhi = kTwoPi * rndArray[2];
    randDirX = randSinTheta * cos(randPhi);
    randDirY = randSinTheta * sin(randPhi);
  }

  // if tehere was any energy deposit add it to parent track doposited energy
  tracks.fEdepV[iintrack] += kerma;

  // loop over the secondaries
  for (int i = 0; i < nSecPart; ++i) {
    if (pid[i] >= TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
      int idummy = pid[i] - 1000000000;
      int Z = idummy / 10000.;
      int A = (idummy - Z * 10000) / 10.;
      double secMass = TPartIndex::I()->GetAprxNuclearMass(Z, A);
      double px = mom[3 * i];
      double py = mom[3 * i + 1];
      double pz = mom[3 * i + 2];
      double secPtot2 = px * px + py * py + pz * pz; // total P^2 [GeV^2]
      tracks.fEdepV[iintrack] += sqrt(secPtot2 + secMass * secMass) - secMass;
      continue;
    }

    int secPDG = TPartIndex::I()->PDG(pid[i]); // Geant V particle code -> particle PGD code
#ifdef USE_VECGEOM_NAVIGATOR
    const Particle *const &secPartPDG = &Particle::GetParticle(secPDG);
#else
    TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
#endif
    double secMass = secPartPDG->Mass();
    double px = mom[3 * i];
    double py = mom[3 * i + 1];
    double pz = mom[3 * i + 2];
    double secPtot2 = px * px + py * py + pz * pz;             // total P^2 [GeV^2]
    double secPtot = sqrt(secPtot2);                     // total P [GeV]
    double secEtot = sqrt(secPtot2 + secMass * secMass); // total energy in [GeV]
    double secEkin = secEtot - secMass;                        // kinetic energy in [GeV]
    // Ekin of the i-th secondary is higher than the threshold
    if (secEkin > energyLimit) { // insert secondary into tracks_v
      GeantTrack &gTrack = td->GetTrack();
      // set the new track properties
      gTrack.fEvent = tracks.fEventV[iintrack];
      gTrack.fEvslot = tracks.fEvslotV[iintrack];
      //       gTrack.fParticle = nTotSecPart;          //index of this particle
      gTrack.fPDG = secPDG;    // PDG code of this particle
      gTrack.fGVcode = pid[i]; // GV index of this particle
      gTrack.fEindex = 0;
      gTrack.fCharge = secPartPDG->Charge() / 3.; // charge of this particle
      gTrack.fProcess = 0;
      gTrack.fVindex = tracks.fVindexV[iintrack]; // volume index
      gTrack.fNsteps = 0;
      //       gTrack.fSpecies  = 0;
      gTrack.fStatus = kNew;                  // status of this particle
      gTrack.fMass = secMass;                 // mass of this particle
      gTrack.fXpos = tracks.fXposV[iintrack]; // rx of this particle (same as parent)
      gTrack.fYpos = tracks.fYposV[iintrack]; // ry of this particle (same as parent)
      gTrack.fZpos = tracks.fZposV[iintrack]; // rz of this particle (same as parent)
      gTrack.fXdir = px / secPtot;            // dirx of this particle (before transform.)
      gTrack.fYdir = py / secPtot;            // diry of this particle before transform.)
      gTrack.fZdir = pz / secPtot;            // dirz of this particle before transform.)
      gTrack.fP = secPtot;                    // momentum of this particle
      gTrack.fE = secEtot;                    // total E of this particle
      gTrack.fTime = tracks.fTimeV[iintrack]; // global time for this particle
      gTrack.fEdep = 0.;
      gTrack.fPstep = 0.;
      gTrack.fStep = 0.;
      gTrack.fSnext = 0.;
      gTrack.fSafety = tracks.fSafetyV[iintrack];
      gTrack.fFrombdr = tracks.fFrombdrV[iintrack];
      gTrack.fPending = kFALSE;
      *gTrack.fPath = *tracks.fPathV[iintrack];
      *gTrack.fNextpath = *tracks.fPathV[iintrack];

      // rotate at-rest secondary by a common random theta and random phi
      RotateNewTrack(randDirX, randDirY, randDirZ, gTrack);

      gPropagator->AddTrack(gTrack);
      tracks.AddTrack(gTrack);

      ++nTotSecPart; // increase # of secondaries in tracks_v
    } else {
      // add the Ekin of this secondary to the energy depositon
      tracks.fEdepV[iintrack] += secEkin;
      // check if it is a stopped particle and call at-rest sampling if necessary
      if (secEkin <= 0.0 && HasRestProcess(pid[i]))
        GetRestFinStates(pid[i], mxs, energyLimit, tracks, iintrack, nTotSecPart, td); // RECURSION
    }
  } // end loop over the secondaries
}

//______________________________________________________________________________
void TTabPhysMgr::SampleDecayInFlight(int partindex, TMXsec *mxs, double energyLimit, GeantTrack_v &tracks,
                                      int iintrack, int &nTotSecPart, GeantTaskData *td) {
  int nSecPart = 0;     // number of secondary particles per reaction
  const int *pid = 0;   // GeantV particle codes [nSecPart]
  const float *mom = 0; // momentum vectors the secondaries [3*nSecPart]
  char isSurv = 0;      // is the primary survived the interaction

  isSurv = fDecay->SampleDecay(partindex, nSecPart, pid, mom);
  // isSurv should always be FALSE here because primary was stopped
  if (isSurv)
    std::cout << "\n---       A particle survived its decay!!!       ---\n"
              << "----    In TTabPhysMgr::SampleFinalStateAtRest     ---\n" << std::endl;

  if (nSecPart) {
    // Go for the secondaries
    double beta = tracks.fPV[iintrack] / tracks.fEV[iintrack];
    double bx = tracks.fXdirV[iintrack] * beta;
    double by = tracks.fYdirV[iintrack] * beta;
    double bz = tracks.fZdirV[iintrack] * beta;
    double b2 = bx * bx + by * by + bz * bz; // it is beta*beta
    double gam = 1.0 / sqrt(1.0 - b2);
    double gam2 = b2 > 0.0 ? (gam - 1.0) / b2 : 0.0;

    for (int isec = 0; isec < nSecPart; ++isec) {
      if (pid[isec] >= TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
        int idummy = pid[isec] - 1000000000;
        int Z = idummy / 10000.;
        int A = (idummy - Z * 10000) / 10.;
        double secMass = TPartIndex::I()->GetAprxNuclearMass(Z, A);
        double px = mom[3 * isec];
        double py = mom[3 * isec + 1];
        double pz = mom[3 * isec + 2];
        double secPtot2 = px * px + py * py + pz * pz; // total P^2 [GeV^2]
        tracks.fEdepV[iintrack] += sqrt(secPtot2 + secMass * secMass) - secMass;
        continue;
      }

      int secPDG = TPartIndex::I()->PDG(pid[isec]); // GV part.code -> PGD code
#ifdef USE_VECGEOM_NAVIGATOR
      const Particle *const &secPartPDG = &Particle::GetParticle(secPDG);
#else
      TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
#endif
      double secMass = secPartPDG->Mass(); // mass [GeV]
      double px = mom[3 * isec];
      double py = mom[3 * isec + 1];
      double pz = mom[3 * isec + 2];
      double secP2 = px * px + py * py + pz * pz;             // total P^2 [GeV^2]
      double secEtot = sqrt(secP2 + secMass * secMass); // total E [GeV]
      // double secEkin  = secEtot - secMass; //kinetic energy in [GeV]

      double bp = bx * px + by * py + bz * pz;
      px = px + gam2 * bp * bx + gam * bx * secEtot;
      py = py + gam2 * bp * by + gam * by * secEtot;
      pz = pz + gam2 * bp * bz + gam * bz * secEtot;
      secEtot = gam * (secEtot + bp);

      double secPtot = sqrt((secEtot - secMass) * (secEtot + secMass));
      double secEkin = secEtot - secMass;
      if (secEkin > energyLimit) { // insert secondary into tracks_v
        GeantTrack &gTrack = td->GetTrack();
        // set the new track properties
        gTrack.fEvent = tracks.fEventV[iintrack];
        gTrack.fEvslot = tracks.fEvslotV[iintrack];
        //         gTrack.fParticle = nTotSecPart;          //index of this particle
        gTrack.fPDG = secPDG;       // PDG code of this particle
        gTrack.fGVcode = pid[isec]; // GV index of this particle
        gTrack.fEindex = 0;
        gTrack.fCharge = secPartPDG->Charge() / 3.; // charge of this particle
        gTrack.fProcess = -1;
        gTrack.fVindex = tracks.fVindexV[iintrack];
        gTrack.fNsteps = 0;
        //         gTrack.fSpecies  = 0;
        gTrack.fStatus = kNew;                  // status of this particle
        gTrack.fMass = secMass;                 // mass of this particle
        gTrack.fXpos = tracks.fXposV[iintrack]; // rx of this particle (same as parent)
        gTrack.fYpos = tracks.fYposV[iintrack]; // ry of this particle (same as parent)
        gTrack.fZpos = tracks.fZposV[iintrack]; // rz of this particle (same as parent)
        gTrack.fXdir = px / secPtot;            // dirx of this particle (before transform.)
        gTrack.fYdir = py / secPtot;            // diry of this particle before transform.)
        gTrack.fZdir = pz / secPtot;            // dirz of this particle before transform.)
        gTrack.fP = secPtot;                    // momentum of this particle
        gTrack.fE = secEtot;                    // total E of this particle
        gTrack.fTime = tracks.fTimeV[iintrack]; // global time for this track
        gTrack.fEdep = 0.;
        gTrack.fPstep = 0.;
        gTrack.fStep = 0.;
        gTrack.fSnext = 0.;
        gTrack.fSafety = tracks.fSafetyV[iintrack];
        gTrack.fFrombdr = tracks.fFrombdrV[iintrack];
        gTrack.fPending = kFALSE;
        *gTrack.fPath = *tracks.fPathV[iintrack];
        *gTrack.fNextpath = *tracks.fPathV[iintrack];

        gPropagator->AddTrack(gTrack);
        tracks.AddTrack(gTrack);

        ++nTotSecPart; // increase # of secondaries in tracks_v
      } else {
        // add the Ekin of this secondary to the energy depositon
        tracks.fEdepV[iintrack] += secEkin;
        // check if it is a stopped particle and call at-rest sampling if necessary
        if (secEkin <= 0.0 && HasRestProcess(pid[isec]))
          GetRestFinStates(pid[isec], mxs, energyLimit, tracks, iintrack, nTotSecPart, td);
      }
    } // end loop over secondaries
  }   // end if has secondaries
}

//_____________________________________________________________________________
// FOR A SINGLE GeantTrack
// (oldXdir, oldYdir, oldZdir) is the direction vector of parent track in lab.
// frame; direction vector of the current track, measured from local Z is
// already updated in GeantTrack track; here we rotate it to lab. frame
void TTabPhysMgr::RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack &track) {
  const double one = 1.0;
  const double zero = 0.0;
  const double amin = 1.0e-10;
  const double one5 = 1.5;
  const double half = 0.5;

  double cosTheta0 = oldZdir;
  double sinTheta0 = sqrt(oldXdir * oldXdir + oldYdir * oldYdir);
  double cosPhi0;
  double sinPhi0;

  if (sinTheta0 > amin) {
    cosPhi0 = oldXdir / sinTheta0;
    sinPhi0 = oldYdir / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = track.fXdir;
  double h1 = sinTheta0 * track.fZdir + cosTheta0 * h0;
  double h2 = track.fYdir;

  track.fXdir = h1 * cosPhi0 - h2 * sinPhi0;
  track.fYdir = h1 * sinPhi0 + h2 * cosPhi0;
  track.fZdir = track.fZdir * cosTheta0 - h0 * sinTheta0;

  // renormalization: -use 1-th order Taylor aprx. of 1/sqrt(x) around 1.0
  // that should be almost exact since the vector almost normalized!
  double delta = one5 - half * (track.fXdir * track.fXdir + track.fYdir * track.fYdir + track.fZdir * track.fZdir);
  track.fXdir *= delta;
  track.fYdir *= delta;
  track.fZdir *= delta;
}

//______________________________________________________________________________
// FOR THE itrack-th element of a GeantTrack_v
// (oldXdir, oldYdir, oldZdir) is the direction vector of parent track in lab.
// frame; direction vector of the current track, measured from local Z is
// already updated in GeantTrack track; here we rotate it to lab. frame
GEANT_CUDA_DEVICE_CODE
void TTabPhysMgr::RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack_v &tracks,
                                 int itrack) {
  const double one = 1.0;
  const double zero = 0.0;
  const double amin = 1.0e-10;
  const double one5 = 1.5;
  const double half = 0.5;

  double cosTheta0 = oldZdir;
  double sinTheta0 = sqrt(oldXdir * oldXdir + oldYdir * oldYdir);
  double cosPhi0;
  double sinPhi0;

  if (sinTheta0 > amin) {
    cosPhi0 = oldXdir / sinTheta0;
    sinPhi0 = oldYdir / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = tracks.fXdirV[itrack];
  double h1 = sinTheta0 * tracks.fZdirV[itrack] + cosTheta0 * h0;
  double h2 = tracks.fYdirV[itrack];

  tracks.fXdirV[itrack] = h1 * cosPhi0 - h2 * sinPhi0;
  tracks.fYdirV[itrack] = h1 * sinPhi0 + h2 * cosPhi0;
  tracks.fZdirV[itrack] = tracks.fZdirV[itrack] * cosTheta0 - h0 * sinTheta0;

  // renormalization: -use 1-th order Taylor aprx. of 1/sqrt(x) around 1.0
  // that should be almost exact since the vector almost normalized!
  double delta =
      one5 -
      half * (tracks.fXdirV[itrack] * tracks.fXdirV[itrack] + tracks.fYdirV[itrack] * tracks.fYdirV[itrack] +
              tracks.fZdirV[itrack] * tracks.fZdirV[itrack]);
  tracks.fXdirV[itrack] *= delta;
  tracks.fYdirV[itrack] *= delta;
  tracks.fZdirV[itrack] *= delta;
}

//______________________________________________________________________________
// FOR A SINGLE GeantTrack
// GeantTrack track contains the original direction in lab frame; theta and
// phi are the scattering angles measured form the particle local Z
void TTabPhysMgr::RotateTrack(GeantTrack &track, double theta, double phi) {
  const double one = 1.0;
  const double zero = 0.0;
  const double amin = 1.0e-10;
  const double one5 = 1.5;
  const double half = 0.5;

  double cosTheta0 = track.fZdir;
  double sinTheta0 = sqrt(track.fXdir * track.fXdir + track.fYdir * track.fYdir);
  double cosPhi0;
  double sinPhi0;
  double cosTheta = cos(theta);
  double sinTheta = sin(theta);

  if (sinTheta0 > amin) {
    cosPhi0 = track.fXdir / sinTheta0;
    sinPhi0 = track.fYdir / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = sinTheta * cos(phi);
  double h1 = sinTheta0 * cosTheta + cosTheta0 * h0;
  double h2 = sinTheta * sin(phi);

  track.fXdir = h1 * cosPhi0 - h2 * sinPhi0;
  track.fYdir = h1 * sinPhi0 + h2 * cosPhi0;
  track.fZdir = cosTheta * cosTheta0 - h0 * sinTheta0;

  // renormalization: -ensure normality to avoid accumulated numerical errors
  //    due to sequential calls of rotation; avoid 1/sqrt(x) computation by
  //    using the 1-th order Taylor aprx. around 1.0 that should be almost
  //    exact since the vector almost normalized!
  double delta = one5 - half * (track.fXdir * track.fXdir + track.fYdir * track.fYdir + track.fZdir * track.fZdir);
  track.fXdir *= delta;
  track.fYdir *= delta;
  track.fZdir *= delta;
}

//______________________________________________________________________________
// FOR THE itrack-th element of a GeantTrack_v
// GeantTrack_v contains the original direction in lab frame; theta and
// phi are the scattering angles measured form the particle local Z
GEANT_CUDA_DEVICE_CODE
void TTabPhysMgr::RotateTrack(GeantTrack_v &tracks, int itrack, double theta, double phi) {
  const double one = 1.0;
  const double zero = 0.0;
  const double amin = 1.0e-10;
  const double one5 = 1.5;
  const double half = 0.5;

  double cosTheta0 = tracks.fZdirV[itrack];
  double sinTheta0 =
      sqrt(tracks.fXdirV[itrack] * tracks.fXdirV[itrack] + tracks.fYdirV[itrack] * tracks.fYdirV[itrack]);
  double cosPhi0;
  double sinPhi0;
  double cosTheta = cos(theta);
  double sinTheta = sin(theta);

  if (sinTheta0 > amin) {
    cosPhi0 = tracks.fXdirV[itrack] / sinTheta0;
    sinPhi0 = tracks.fYdirV[itrack] / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = sinTheta * cos(phi);
  double h1 = sinTheta0 * cosTheta + cosTheta0 * h0;
  double h2 = sinTheta * sin(phi);

  tracks.fXdirV[itrack] = h1 * cosPhi0 - h2 * sinPhi0;
  tracks.fYdirV[itrack] = h1 * sinPhi0 + h2 * cosPhi0;
  tracks.fZdirV[itrack] = cosTheta * cosTheta0 - h0 * sinTheta0;

  // renormalization: -ensure normality to avoid accumulated numerical errors
  //    due to sequential calls of rotation; avoid 1/sqrt(x) computation by
  //    using the 1-th order Taylor aprx. around 1.0 that should be almost
  //    exact since the vector almost normalized!
  double delta =
      one5 -
      half * (tracks.fXdirV[itrack] * tracks.fXdirV[itrack] + tracks.fYdirV[itrack] * tracks.fYdirV[itrack] +
              tracks.fZdirV[itrack] * tracks.fZdirV[itrack]);
  tracks.fXdirV[itrack] *= delta;
  tracks.fYdirV[itrack] *= delta;
  tracks.fZdirV[itrack] *= delta;
}

//______________________________________________________________________________
char *TTabPhysMgr::GetVersion() {
  char *ver = new char[512];
  sprintf(ver, "%d.%d.%d", VersionMajor(), VersionMinor(), VersionSub());
  return ver;
}

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
bool TTabPhysMgr::HasRestProcess(int gvindex) {
  return fDecay->HasDecay(gvindex) || fHasNCaptureAtRest[gvindex] || (gvindex == TPartIndex::I()->GetSpecGVIndex(1));
}
