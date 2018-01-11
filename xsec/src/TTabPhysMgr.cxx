#include "TTabPhysMgr.h"

#include <VecCore/VecCore>

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#include "base/RNG.h"
#include "Types.h"
#include "SystemOfUnits.h"
#include "Material.h"
#include "MaterialProperties.h"
#include "Element.h"
#else
#include "TGeoManager.h"
#include "TGeoBranchArray.h"
#include "TGeoExtension.h"
#include "TList.h"
#endif

#include "GeantTrackVec.h"

#include "globals.h"
#include "GeantPropagator.h"
#include "GeantTaskData.h"

#ifdef USE_ROOT
#include "TRandom.h"
#include "TFile.h"
#include "TBits.h"
#include "TError.h"
#include "TSystem.h"
#endif

#include "TPartIndex.h"
#include "TEXsec.h"
#include "TMXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"

#include <iostream>
#include <time.h>

#include "Geant/Math.h"
using vecgeom::kPi;
using vecgeom::kTwoPi;
#ifndef VECCORE_CUDA
#include "base/MessageLogger.h"
#endif
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

TTabPhysMgr *TTabPhysMgr::fgInstance = 0;
using namespace Geant;

//______________________________________________________________________________
TTabPhysMgr *TTabPhysMgr::Instance(const char *xsecfilename, const char *finalsfilename) {
  // Access to instance of TTabPhysMgr
#ifndef VECCORE_CUDA
  if (fgInstance)
    return fgInstance;
  if (!(xsecfilename && finalsfilename)) {
    log_error(std::cout, "TTabPhysMgr::Instance", "Create TTabPhysMgr instance providing xsec files\n");
    return 0;
  }
  fgInstance = new TTabPhysMgr(xsecfilename, finalsfilename);
  return fgInstance;
#else
  (void)xsecfilename;
  (void)finalsfilename;
  // ... missing return ...
  return nullptr;
#endif
}

//______________________________________________________________________________
TTabPhysMgr::~TTabPhysMgr() {
  // Destructor
#ifndef VECCORE_CUDA
  delete[] fMatXsec;
  delete[] fElemXsec;
  delete[] fElemFstate;
  delete fDecay;
  delete fHasNCaptureAtRest;
  fgInstance = 0;
//  Particle_t::CreateParticles();
#endif
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
#ifndef VECCORE_CUDA
  fgInstance = this;
  //Particle_t::CreateParticles();
#ifdef USE_ROOT
  clock_t t = clock();
// Load elements from geometry, however in most cases it should already be done
#ifdef USE_VECGEOM_NAVIGATOR
  //td->fPropagator->LoadVecGeomGeometry();
  const geantphysics::Vector_t<Material_t*> matlist = geantphysics::Material::GetTheMaterialTable();
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
  fDecay = TEFstate::GetDecayTable();
  if(fDecay == nullptr) fDecay = (TPDecay *)fstate->Get("DecayTable");

#ifdef USE_VECGEOM_NAVIGATOR
  printf("#materials:= %lu \n", matlist.size());
#else
  // INFO: print number of materials in the current GeoManager
  printf("#materials:= %d \n", matlist->GetSize());
#endif

  // First loop on all materials to mark used elements
  TBits elements(NELEM);
#ifdef USE_VECGEOM_NAVIGATOR
  for (unsigned int i = 0; i < matlist.size(); ++i) {
    mat = matlist[i];
    double zeff = mat->GetMaterialProperties()->GetEffectiveZ();
    std::cout << mat->GetName() << "used " << mat->IsUsed() << " Zeff = "<< zeff << std::endl;
    if (!mat->IsUsed() || zeff < 1.)
      continue;
    fNmaterials++;
    int nelem = mat->GetNumberOfElements();
    // Check if we are on the safe side; should exit otherwise
    if (nelem > MAXNELEMENTS) {
      Fatal("TTabPhysMgr", "Number of elements in %s is %d > TTabPhysMgr::MAXNELEMENTS=%d\n", mat->GetName().c_str(), nelem,
            MAXNELEMENTS);
    }
    const geantphysics::Vector_t<geantphysics::Element*> elemVect =  mat->GetElementVector();
    for (int iel=0; iel<nelem; ++iel) {
      double zd = elemVect[iel]->GetZ();
      if (zd < 1 || zd > NELEM) {
        Fatal("TTabPhysMgr", "In material %s found element with z=%d > NELEM=%d", mat->GetName().c_str(), (int)zd, NELEM);
      }
      elements.SetBitNumber(zd);
    }
  }
#else
  while ((mat = (Material_t *)next())) {
    std::cout << mat->GetName() << "used " << mat->IsUsed() << " Zeff = "<< mat->GetZ() << std::endl;
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
#endif

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

#ifdef USE_VECGEOM_NAVIGATOR
  int *z = new int[MAXNELEMENTS];
  int *a = new int[MAXNELEMENTS];
  float *w = new float[MAXNELEMENTS];
  fNmaterials = 0;
  for (unsigned int i = 0; i < matlist.size(); ++i) {
    mat = matlist[i];
    if (!mat->IsUsed())
      continue;
    const geantphysics::Vector_t<geantphysics::Element*> theElemVect = mat->GetElementVector();
    const double *massFractionVect = mat->GetMassFractionVector();
    int nelem = mat->GetNumberOfElements();
//    std::cout<< " ==== material = " << mat->GetName() << std::endl;
    for (int iel = 0; iel < nelem; ++iel) {
       z[iel] = theElemVect[iel]->GetZ();
       a[iel] = theElemVect[iel]->GetA()/(geant::g/geant::mole);
       w[iel] = massFractionVect[iel];
//       std::cout<< "iel = "<< iel <<" z = " << z[iel] << " a = " << a[iel] << " w = " << w[iel] << " density = " << mat->GetDensity()/(geant::g/geant::cm3) << std::endl;
    }
    if (nelem == 0) {
      std::cout<<mat<<std::endl;
      Fatal("TTabPhysMgr", "The material (%s) seems to have no elements", mat->GetName().c_str());
    }
    // Construct the TMXsec object that corresponds to the current material
    TMXsec *mxs = new TMXsec(mat->GetName().c_str(), mat->GetName().c_str(), z, a, w, nelem, mat->GetDensity()/(geant::g/geant::cm3), true, fDecay);
    fMatXsec[fNmaterials++] = mxs;
    mat->SetXsecPtr(static_cast<void *>(mxs));
  } // End of while
  delete[] z;
  delete[] a;
  delete[] w;
#else
  int *z = new int[MAXNELEMENTS];
  int *a = new int[MAXNELEMENTS];
  float *w = new float[MAXNELEMENTS];
  fNmaterials = 0;
  next.Reset();
  while ((mat = (Material_t *)next())) {
    if (!mat->IsUsed())
      continue;
    int nelem = mat->GetNelements();
    // loop over the elements of the current material in order to obtain the
    // z, a, w, arrays of the elements of this material
    double ad;
    double zd;
    double wd;
//    std::cout<< " ==== material = " << mat->GetName() << std::endl;
    for (int iel = 0; iel < nelem; ++iel) {
      mat->GetElementProp(ad, zd, wd, iel);
      a[iel] = ad;
      z[iel] = zd;
      w[iel] = wd;
//      std::cout<< "iel = "<< iel <<" z = " << z[iel] << " a = " << a[iel] << " w = " << w[iel] << " density = " << mat->GetDensity() << std::endl;
    }
    if (nelem == 0) {
      mat->Dump();
      Fatal("TTabPhysMgr", "The material (%s) seems to have no elements", mat->GetName());
    }
    // Construct the TMXsec object that corresponds to the current material
    TMXsec *mxs = new TMXsec(mat->GetName(), mat->GetName(), z, a, w, nelem, mat->GetDensity(), true, fDecay);
    fMatXsec[fNmaterials++] = mxs;
    mat->SetFWExtension(new TGeoRCExtension(new TOMXsec(mxs)));
  } // End of while
  delete[] z;
  delete[] a;
  delete[] w;
#endif

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
  t = clock() - t;
  printf("Memory taken by xsec and states: %ld [MB] loaded in: %g [sec]\n", mem, ((float)t) / CLOCKS_PER_SEC);
#endif
#endif
  (void)xsecfilename;
  (void)finalsfilename;
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

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::TransformLF(int /*indref*/, TrackVec_t & /*tracks*/, int /*nproducts*/, int /*indprod*/,
                              TrackVec_t & /*output*/) {
  // Transform tracks taken from the final state from the local frame to the lab
  // frame (LF). Not clear what parameters to add yet.
  // Input: reference track (mother) described as vector container + index of ref track
  // Input: number of tracks in the final state, start index and vector container
  // Output: roto-boosted tracks in the output vector
}

// NOT ACTIVE NOW
//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::ApplyMsc(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Compute MSC angle at the beginning of the step and apply it to the vector
  // of tracks.
  // Input: material index, number of tracks in the tracks vector to be used
  // Output: fXdirV, fYdirV, fZdirV modified in the track container for ntracks

  TMXsec *mxs = 0;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
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

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  double *rndArray = td->GetDblArray(ntracks);
#ifdef USE_VECGEOM_NAVIGATOR
  td->fRndm->uniform_array(ntracks, rndArray);
#elif USE_ROOT
  td->fRndm->RndmArray(ntracks, rndArray);
#endif
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
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::ApplyMsc(Material_t *mat, TrackVec_t &tracks, GeantTaskData *td) {
  // Compute MSC angle at the beginning of the step and apply it to the vector
  // of tracks.
  // Input: material index, number of tracks in the tracks vector to be used
  // Output: fXdirV, fYdirV, fZdirV modified in the track container for ntracks

  TMXsec *mxs = 0;
  int ntracks = tracks.size();
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
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

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  double *rndArray = td->GetDblArray(ntracks);
#ifdef USE_VECGEOM_NAVIGATOR
  td->fRndm->uniform_array(ntracks, rndArray);
#elif USE_ROOT
  td->fRndm->RndmArray(ntracks, rndArray);
#endif
#else
  double *rndArray = 0; // NOTE: we need to get it from somewhere ....
  VECGEOM_NAMESPACE::RNG::Instance().uniform_array(ntracks, rndArray, 0., 1.);
#endif

  //   double dir[3] = {0.,0.,0.};
  if (mxs) {
    for (int i = 0; i < ntracks; ++i) {
      msTheta = mxs->MS(tracks[i]->GVcode(), tracks[i]->T());
      msPhi = 2. * kPi * rndArray[i];
      RotateTrack(*tracks[i], msTheta, msPhi);
    }
    return;
  }
  // Mixed tracks in different volumes
  for (int i = 0; i < ntracks; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)tracks[i]->GetMaterial()->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)tracks[i]->GetMaterial()->GetFWExtension())->GetUserObject())->MXsec();
#endif
    msTheta = mxs->MS(tracks[i]->GVcode(), tracks[i]->T());
    msPhi = 2. * kPi * rndArray[i];
    /*
          if (icnt<100 && mat->GetZ()>10) {
             Printf("theta=%g  phi=%g", msTheta*vecgeom::Materialh::RadToDeg(), msPhi*vecgeom::Materialh::RadToDeg());
             dir[0] = tracks[i]->fXdir;
             dir[1] = tracks[i]->fYdir;
             dir[2] = tracks[i]->fZdir;
          }
    */
    RotateTrack(*tracks[i], msTheta, msPhi);
    /*
          if (icnt<100 && mat->GetZ()>10) {
             icnt++;
             double dot = dir[0]*tracks[i]->fXdir + dir[1]*tracks[i]->fYdir +dir[2]*tracks[i]->fZdir;
             double angle = vecgeom::Materialh::ACos(dot)*vecgeom::Materialh::RadToDeg();
             Printf("new angle=%g   delta=%g", angle,
       vecgeom::Materialh::Abs(angle-msTheta*vecgeom::Materialh::RadToDeg()));
          }
    */
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TTabPhysMgr::Eloss(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Apply energy loss for the input material for ntracks in the vector of
  // tracks. Output: modified tracks.fEV array

  TMXsec *mxs = 0;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
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
  double energyLimit = td->fPropagator->fConfig->fEmin;
  if (mxs) {
    mxs->Eloss(ntracks, tracks,td);
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
    mxs->ElossSingle(i, tracks,td);
    // call atRest sampling for tracks that have been stopped by Eloss and has at-rest
    if (tracks.fProcessV[i] == -2 && HasRestProcess(tracks.fGVcodeV[i]))
      GetRestFinStates(tracks.fGVcodeV[i], mxs, energyLimit, tracks, i, nTotSecPart, td);
  }

  return nTotSecPart;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TTabPhysMgr::Eloss(GeantTrack *track, TrackVec_t &output, GeantTaskData *td) {
  // Apply energy loss for the input material for ntracks in the vector of
  // tracks. Output: modified tracks.fEV array

  TMXsec *mxs = 0;
  int nTotSecPart = 0; // total number of new tracks
  double energyLimit = td->fPropagator->fConfig->fEmin;
#ifdef USE_VECGEOM_NAVIGATOR
  mxs = (TMXsec *)track->GetMaterial()->GetXsecPtr();
#else
  mxs = ((TOMXsec *)((TGeoRCExtension *)track->GetMaterial()->GetFWExtension())->GetUserObject())->MXsec();
#endif
  mxs->Eloss(*track, td);
  // call atRest sampling for tracks that have been stopped by Eloss and has at-rest
  if (track->Process() == -2 && HasRestProcess(track->GVcode()))
    GetRestFinStates(track->GVcode(), mxs, energyLimit, track, nTotSecPart, output, td);

  return nTotSecPart;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TTabPhysMgr::Eloss(TrackVec_t &tracks, TrackVec_t &output, GeantTaskData *td) {
  // Apply energy loss for the input material for ntracks in the vector of
  // tracks. Output: modified tracks.fEV array

  TMXsec *mxs = 0;
  int ntracks = tracks.size();
  int nTotSecPart = 0; // total number of new tracks
  double energyLimit = td->fPropagator->fConfig->fEmin;
  // Mixed tracks in different volumes
  for (int i = 0; i < ntracks; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)tracks[i]->GetMaterial()->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)tracks[i]->GetMaterial()->GetFWExtension())->GetUserObject())->MXsec();
#endif
    mxs->Eloss(*tracks[i], td);
    // call atRest sampling for tracks that have been stopped by Eloss and has at-rest
    if (tracks[i]->Process() == -2 && HasRestProcess(tracks[i]->GVcode()))
      GetRestFinStates(tracks[i]->GVcode(), mxs, energyLimit, tracks[i], nTotSecPart, output, td);
  }

  return nTotSecPart;
}

//______________________________________________________________________________
void TTabPhysMgr::ProposeStep(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Sample free flight/proposed step for the firts ntracks tracks and store them
  // in tracks.fPstepV

  TMXsec *mxs = 0;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (mat) {
#ifdef USE_VECGEOM_NAVIGATOR
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

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::ProposeStep(GeantTrack *track, GeantTaskData *td) {
  // Sample free flight/proposed step for the firts ntracks tracks and store them
  // in tracks.fPstepV

  TMXsec *mxs = 0;
#ifdef USE_VECGEOM_NAVIGATOR
  mxs = (TMXsec *)track->GetMaterial()->GetXsecPtr();
#else
  mxs = ((TOMXsec *)((TGeoRCExtension *)track->GetMaterial()->GetFWExtension())->GetUserObject())->MXsec();
#endif
  mxs->ProposeStep(*track, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::ProposeStep(TrackVec_t &tracks, GeantTaskData *td) {
  // Sample free flight/proposed step for the firts ntracks tracks and store them
  // in tracks.fPstepV

  TMXsec *mxs = 0;
  int ntracks = tracks.size();
  for (int i = 0; i < ntracks; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)tracks[i]->GetMaterial()->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)tracks[i]->GetMaterial()->GetFWExtension())->GetUserObject())->MXsec();
#endif
    mxs->ProposeStep(*tracks[i], td);
  }
}

// Implemented in a different way
//______________________________________________________________________________
int TTabPhysMgr::SampleDecay(int /*ntracks*/, GeantTrack_v & /*tracksin*/, GeantTrack_v & /*tracksout*/) {
  // Sample decay for the tracks in the input vector and push the resulting tracks in
  // the output vector. Change status of decayed tracks. Returns number of new tracks.
  return 0;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TTabPhysMgr::SampleDecay(TrackVec_t &/*tracksin*/, TrackVec_t &/*tracksout*/) {
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
    mat = Material_t::GetTheMaterialTable()[imat];
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

// sampling: target atom and type of the interaction for each primary tracks
//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::SampleTypeOfInteractions(GeantTrack *track, GeantTaskData *td) {
  TMXsec *mxs = 0;
#ifdef USE_VECGEOM_NAVIGATOR
  mxs = (TMXsec *)track->GetMaterial()->GetXsecPtr();
#else
  mxs = ((TOMXsec *)((TGeoRCExtension *)track->GetMaterial()->GetFWExtension())->GetUserObject())->MXsec();
#endif
  mxs->SampleInt(*track, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::SampleTypeOfInteractions(TrackVec_t &tracks, GeantTaskData *td) {
  TMXsec *mxs = 0;
  int ntracks = tracks.size();
  for (int i = 0; i < ntracks; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
    mxs = (TMXsec *)tracks[i]->GetMaterial()->GetXsecPtr();
#else
    mxs = ((TOMXsec *)((TGeoRCExtension *)tracks[i]->GetMaterial()->GetFWExtension())->GetUserObject())->MXsec();
#endif
    mxs->SampleInt(*tracks[i], td);
  }
}

//______________________________________________________________________________
int TTabPhysMgr::SampleFinalStates(int imat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  GeantPropagator *propagator = td->fPropagator;
  double energyLimit = propagator->fConfig->fEmin;

  Material_t *mat = 0;
  TMXsec *mxs = 0;
  if (imat >= 0) {
#ifdef USE_VECGEOM_NAVIGATOR
    mat = Material_t::GetTheMaterialTable()[imat];
    mxs = (TMXsec *)mat->GetXsecPtr();
#else
    mat = (Material_t *)fGeom->GetListOfMaterials()->At(imat);
    mxs = ((TOMXsec *)((TGeoRCExtension *)mat->GetFWExtension())->GetUserObject())->MXsec();
#endif
  }

  // tid-based rng
  double *rndArray = td->GetDblArray(2*ntracks);
#ifdef USE_VECGEOM_NAVIGATOR
  td->fRndm->uniform_array(2 * ntracks, rndArray);
#elif USE_ROOT
  td->fRndm->RndmArray(2 * ntracks, rndArray);
#endif

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

    // first check the results of interaction sampling:
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
          vecCore::math::Sqrt(curPrimEkin * (curPrimEkin + 2.0 * primMass) / (orgPrimEkin * (orgPrimEkin + 2.0 * primMass)));
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
        postEkinOfParimary = vecCore::math::Sqrt(postPrimP2 + primMass * primMass) - primMass;
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
        postEkinOfParimary = vecCore::math::Sqrt(postPrimP2 + primMass * primMass) - primMass;

        // update primary in tracks
        double secPtot = vecCore::math::Sqrt(postPrimP2);                      // total P [GeV]
        double secEtot = postEkinOfParimary + tracks.fMassV[t]; // total energy in [GeV]
        tracks.fPV[t] = secPtot;                                // momentum of this particle
        tracks.fEV[t] = secEtot;                                // total E of this particle
        tracks.fXdirV[t] = px / secPtot;                        // dirx of this particle (before transform.)
        tracks.fYdirV[t] = py / secPtot;                        // diry of this particle (before transform.)
        tracks.fZdirV[t] = pz / secPtot;                        // dirz of this particle (before transform.)

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
      // j=0 -> including stopped primary as well if isSurv = true;
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
          tracks.fEdepV[t] += vecCore::math::Sqrt(secPtot2 + secMass * secMass) - secMass;
          continue;
        }
        int secPDG = TPartIndex::I()->PDG(pid[i]); // Geant V particle code -> particle PGD code
#ifdef USE_VECGEOM_NAVIGATOR
        const Particle_t *const &secPartPDG = &Particle_t::GetParticle(secPDG);
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
        double secPtot2 = px * px + py * py + pz * pz;       // total P^2 [GeV^2]
        double secEtot = vecCore::math::Sqrt(secPtot2 + secMass * secMass); // total energy in [Ge
        double secEkin = secEtot - secMass;                  // kinetic energy in [GeV]
        // Ekin of the i-th secondary is higher than the threshold
        if (secEkin >= energyLimit) { // insert secondary into OUT tracks_v and rotate
          GeantTrack &track = td->GetNewTrack();
          double secPtot = vecCore::math::Sqrt(secPtot2);                     // total P [GeV]
          double inv_secPtot = 1.0 / secPtot;
          //          GeantTrack track;
          // set the new track properties
          track.SetEvent(tracks.fEventV[t]);
          track.SetEvslot(tracks.fEvslotV[t]);
          track.SetPDG(secPDG);    // PDG code of this particle
          track.SetGVcode(pid[i]); // GV index of this particle
          track.SetCharge(secPartPDG->Charge()); // charge of this particle
#ifndef USE_VECGEOM_NAVIGATOR
	        track.SetCharge(track.Charge()/3.);
#endif
          track.SetStatus(kNew);           // status of this particle
          track.SetMass(secMass);          // mass of this particle
          track.SetPosition(tracks.fXposV[t], tracks.fYposV[t], tracks.fZposV[t]);
          track.SetDirection(px * inv_secPtot, py * inv_secPtot, pz * inv_secPtot);
          track.SetP(secPtot);             // momentum of this particle
          track.SetE(secEtot);             // total E of this particle
          track.SetTime(tracks.fTimeV[t]); // global time
          track.SetSafety(tracks.fSafetyV[t]);
          track.SetBoundary(tracks.fBoundaryV[t]);
          track.SetProcess(tracks.fProcessV[t]); // Record id of creating process -- Was 0
          track.SetPath(tracks.fPathV[t]);
          track.SetNextPath(tracks.fPathV[t]);
          track.SetMother(tracks.fParticleV[t]);

          // Rotate new track to parent track's frame
          RotateNewTrack(oldXdir, oldYdir, oldZdir, track);

          propagator->AddTrack(track);
          tracks.AddTrack(track);
          td->ReleaseTrack(track);

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

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TTabPhysMgr::SampleFinalStates(TrackVec_t &tracks, TrackVec_t &output, GeantTaskData *td) {

  int nTotSecPart = 0; // total number of secondary particles in tracks
  int ntracks = tracks.size();
  for (int t = 0; t < ntracks; ++t)
    nTotSecPart += SampleFinalStates(tracks[t], output, td);
  return nTotSecPart;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TTabPhysMgr::SampleFinalStates(GeantTrack *track, TrackVec_t &output, GeantTaskData *td) {
  GeantPropagator *propagator = td->fPropagator;
  double energyLimit = propagator->fConfig->fEmin;
  TMXsec *mxs = 0;

  // tid-based rng
  double *rndArray = td->fDblArray;
#ifdef USE_VECGEOM_NAVIGATOR
  td->fRndm->uniform_array(2 , rndArray);
#elif USE_ROOT
  td->fRndm->RndmArray(2 , rndArray);
#endif

  int nTotSecPart = 0; // total number of secondary particles in tracks

  // if no interaction was selected for this track (because it doesn't have any)
  if (track->Process() < 0)
    return 0;

  int nSecPart = 0;     // number of secondary particles per reaction
  const int *pid = 0;   // GeantV particle codes [nSecPart]
  const float *mom = 0; // momentum vectors the secondaries [3*nSecPart]
  float ener = 0;       // energy at the fstate (Ekin of primary after the interc.)
  float kerma = 0;      // released energy
  float weight = 0;     // weight of the fstate (just a dummy parameter now)
  char isSurv = 0;      // is the primary survived the interaction
  int ebinindx = -1;    // energy bin index of the selected final state

#ifdef USE_VECGEOM_NAVIGATOR
  mxs = (TMXsec *)track->GetMaterial()->GetXsecPtr();
#else
  mxs = ((TOMXsec *)((TGeoRCExtension *)track->GetMaterial()->GetFWExtension())->GetUserObject())->MXsec();
#endif

  // first check the results of interaction sampling:
  if (track->Process() == 3) {
    // decay : in-flight decay was selected
    // kill the primary tarck
    track->SetStatus(kKilled);
    // sample in-flight decay final state
    SampleDecayInFlight(track->GVcode(), mxs, energyLimit, track, nTotSecPart, output, td);
    track->SetP(0.);
    track->SetE(track->Mass());
    return nTotSecPart;
  }

  // not decay but something else was selected
  double curPrimEkin = track->T();
  isSurv = fElemFstate[track->EIndex()]->SampleReac(track->GVcode(), track->Process(), curPrimEkin, nSecPart,
                                                   weight, kerma, ener, pid, mom, ebinindx, rndArray[0], rndArray[1]);

  // it is the case of: pre-step energy sigma is not zero of this interaction
  //                    but post step is zero-> we don't have final state for
  //                    this interaction at the postStep energy bin-> do nothing
  // This can happen if interaction type is selected based on the pre-step energy
  // and there is some energy loss along the step (or something like this)
  if (isSurv && ener < 0) // let it go further as it is
    return nSecPart;

  // we should correct the kerma as well but we don't have enough information
  track->IncreaseEdep(kerma);

  // if we have secondaries from the current interaction
  if (nSecPart) {
    double oldXdir = track->Dx(); // old X direction of the primary
    double oldYdir = track->Dy(); // old Y direction of the primary
    double oldZdir = track->Dz(); // old Z direction of the primary
    int j = 0;

    // setting the final state correction factor (we scale only the 3-momentums)
    //-get mass of the primary
    double primMass = track->Mass(); // mass [GeV]
    //-compute corFactor = P_current/P_original = Pz_current/Pz_original
    // (normaly a check would be good but not necessary: if(ebinindx<0 -> ...)
    double orgPrimEkin = (TPartIndex::I()->EGrid())[ebinindx];
    double corFactor =
      vecCore::math::Sqrt(curPrimEkin * (curPrimEkin + 2.0 * primMass) / (orgPrimEkin * (orgPrimEkin + 2.0 * primMass)));
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
      postEkinOfParimary = vecCore::math::Sqrt(postPrimP2 + primMass * primMass) - primMass;
    }

    if (postEkinOfParimary > energyLimit) { // survived even after the correction and the E-limit.
      // keep alive
      //        tracks[t]->fStatus = kAlive;
      double px = mom[0];
      double py = mom[1];
      double pz = mom[2];
      px *= corFactor;
      py *= corFactor;
      pz *= corFactor;
      // compute corrected P^2 in [GeV^2]
      double postPrimP2 = px * px + py * py + pz * pz;
      // recompute post-interaction Ekin of the primary with corrected 3-momentum
      postEkinOfParimary = vecCore::math::Sqrt(postPrimP2 + primMass * primMass) - primMass;

      // update primary in tracks
      double secPtot = vecCore::math::Sqrt(postPrimP2);                      // total P [GeV]
      double secEtot = postEkinOfParimary + track->Mass(); // total energy in [GeV]
      track->SetP(secPtot);                                // momentum of this particle
      track->SetE(secEtot);                                // total E of this particle
      track->SetDirection(px / secPtot, py / secPtot, pz / secPtot);

      // Rotate parent track in tracks to original parent track's frame
      RotateNewTrack(oldXdir, oldYdir, oldZdir, *track);
      // primary track is updated
    } else {
      // Primary particle energy is below tracking limit
      //-set status of primary in tracks to kKilled;
      track->SetStatus(kKilled);
      track->IncreaseEdep(postEkinOfParimary);
      track->SetP(0.);
      track->SetE(track->Mass());
      // if the primary is stopped i.e. Ekin <= 0 then call at-rest if it has
      if (isSurv && postEkinOfParimary <= 0.0 && HasRestProcess(track->GVcode()))
        GetRestFinStates(track->GVcode(), mxs, energyLimit, track, nTotSecPart, output, td);
    }

    if (isSurv)
      j = 1;
    // loop over the secondaries and put them into tracks if they good to track:
    // j=0 -> including stopped primary as well if isSurv = true;
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
        track->IncreaseEdep(vecCore::math::Sqrt(secPtot2 + secMass * secMass) - secMass);
        continue;
      }
      int secPDG = TPartIndex::I()->PDG(pid[i]); // Geant V particle code -> particle PGD code
#ifdef USE_VECGEOM_NAVIGATOR
      const Particle_t *const &secPartPDG = &Particle_t::GetParticle(secPDG);
#else
      TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
#endif
      double secMass = secPartPDG->Mass();
      double px = mom[3 * i];
      double py = mom[3 * i + 1];
      double pz = mom[3 * i + 2];
      px *= corFactor;
      py *= corFactor;
      pz *= corFactor;
      double secPtot2 = px * px + py * py + pz * pz;       // total P^2 [GeV^2]
      double secPtot = vecCore::math::Sqrt(secPtot2);                     // total P [GeV]
      double secEtot = vecCore::math::Sqrt(secPtot2 + secMass * secMass); // total energy in [GeV]
      double secEkin = secEtot - secMass;                  // kinetic energy in [GeV]
      // Ekin of the i-th secondary is higher than the threshold
      if (secEkin >= energyLimit) { // insert secondary into OUT tracks_v and rotate
        GeantTrack &track1 = td->GetNewTrack();
        //          GeantTrack track;
        // set the new track properties
        track1.SetEvent(track->Event());
        track1.SetEvslot(track->EventSlot());
        //          track1.fParticle = nTotSecPart;          //index of this particle
        track1.SetPDG(secPDG);    // PDG code of this particle
        track1.SetGVcode(pid[i]); // GV index of this particle
        track1.SetCharge(secPartPDG->Charge()); // charge of this particle
#ifndef USE_VECGEOM_NAVIGATOR
        track1.SetCharge(track1.Charge()/3.);
#endif
        track1.SetGeneration(track->GetGeneration() + 1);
        //          track1.fSpecies  = 0;
        track1.SetStatus(kNew);           // status of this particle
        track1.SetStage(int(kPreStepStage));
        track1.SetMass(secMass);          // mass of this particle
        track1.SetPosition(track->X(), track->Y(), track->Z());
        track1.SetDirection(px / secPtot, py / secPtot, pz / secPtot);
        track1.SetP(secPtot);             // momentum of this particle
        track1.SetE(secEtot);             // total E of this particle
        track1.SetTime(track->Time()); // global time
        track1.SetSafety(track->GetSafety());
        track1.SetBoundary(track->Boundary());
        track1.SetPath(track->Path());
        track1.SetNextPath(track->Path());
	      track1.SetMother(track->Particle());

        // Rotate new track to parent track's frame
        RotateNewTrack(oldXdir, oldYdir, oldZdir, track1);

        propagator->AddTrack(track1);
        output.push_back(&track1);

        ++nTotSecPart;
      } else {                       // {secondary Ekin < energyLimit} -> kill this secondary
        track->IncreaseEdep(secEkin); // add the Ekin of this secondary to the energy depositon
        // is secEkin <=0 then call at-rest process if the sec. particle has any
        if (secEkin <= 0.0 && HasRestProcess(pid[i]))
          GetRestFinStates(pid[i], mxs, energyLimit, track, nTotSecPart, output, td);
      }
    }                             // end loop over the secondaries
  } else {                        // nSecPart = 0 i.e. there is no any secondaries -> primary was killed as well
    track->SetStatus(kKilled); // set status of primary in tracks to kKilled;
    track->SetP(0);
    track->SetE(track->Mass());
  }

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

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TTabPhysMgr::SampleInt(int imat, TrackVec_t &tracks, GeantTaskData *td)
{
}

*/

// Will be called only if the particle has decay or/and nuclear capture at-rest
//______________________________________________________________________________
// will be called recursively if necessary
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::GetRestFinStates(int partindex, TMXsec *mxs, double energyLimit, GeantTrack_v &tracks, int iintrack,
                                   int &nTotSecPart, GeantTaskData *td) {
  // current track should have already been killed before calling
  const double mecc = 0.00051099906; // e- mass c2 in [GeV]
  double rndArray[3];
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
#ifdef USE_VECGEOM_NAVIGATOR
  td->fRndm->uniform_array(3, rndArray);
#elif USE_ROOT
  td->fRndm->RndmArray(3, rndArray);
#endif
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
      double randSinTheta = vecCore::math::Sqrt(1.0 - randDirZ * randDirZ);
      double randPhi = 2.0 * rndArray[1] * kPi;
      double randDirX = randSinTheta * vecCore::math::Cos(randPhi);
      double randDirY = randSinTheta * vecCore::math::Sin(randPhi);

      // need to do it one-by-one
      // 1. gamma
      GeantTrack &track1 = td->GetNewTrack();
      // set the new track properties: 2 gamma with m_{e}*c*c
      track1.SetEvent(tracks.fEventV[iintrack]);
      track1.SetEvslot(tracks.fEvslotV[iintrack]);
      track1.SetPDG(22);                                    // gamma PDG code
      track1.SetGVcode(TPartIndex::I()->GetSpecGVIndex(2)); // gamma GV index
      track1.SetStatus(kNew);                  // status of this particle
      track1.SetStage(int(kPreStepStage));
      track1.SetPosition(tracks.fXposV[iintrack], tracks.fYposV[iintrack], tracks.fZposV[iintrack]);
      track1.SetDirection(randDirX, randDirY, randDirZ);
      track1.SetP(mecc);                       // momentum of this particle
      track1.SetE(mecc);                       // total E of this particle
      track1.SetTime(tracks.fTimeV[iintrack]); // total time of this particle
      track1.SetSafety(tracks.fSafetyV[iintrack]);
      track1.SetBoundary(tracks.fBoundaryV[iintrack]);
      track1.SetPath(tracks.fPathV[iintrack]);
      track1.SetNextPath(tracks.fPathV[iintrack]);

      td->fPropagator->AddTrack(track1);
      tracks.AddTrack(track1);

      // 2. gamma : everything is the same but the direction
      track1.SetDirection(-randDirX, -randDirY, -randDirZ);

      td->fPropagator->AddTrack(track1);
      tracks.AddTrack(track1);
      td->ReleaseTrack(track1);

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

  // isSurv should always be false here because primary was stopped -> just a check
  if (isSurv)
    printf("A stopped particle survived its rest process in TTabPhysMgr::GetRestFinSates!\n");

  // for a random rotation
  if (nSecPart) {
    randDirZ = 1.0 - 2.0 * rndArray[1];
    randSinTheta = vecCore::math::Sqrt((1.0 - randDirZ) * (1.0 + randDirZ));
    randPhi = kTwoPi * rndArray[2];
    randDirX = randSinTheta * vecCore::math::Cos(randPhi);
    randDirY = randSinTheta * vecCore::math::Sin(randPhi);
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
      tracks.fEdepV[iintrack] += vecCore::math::Sqrt(secPtot2 + secMass * secMass) - secMass;
      continue;
    }

    int secPDG = TPartIndex::I()->PDG(pid[i]); // Geant V particle code -> particle PGD code
#ifdef USE_VECGEOM_NAVIGATOR
    const Particle_t *const &secPartPDG = &Particle_t::GetParticle(secPDG);
#else
    TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
#endif
    double secMass = secPartPDG->Mass();
    double px = mom[3 * i];
    double py = mom[3 * i + 1];
    double pz = mom[3 * i + 2];
    double secPtot2 = px * px + py * py + pz * pz;       // total P^2 [GeV^2]
    double secPtot = vecCore::math::Sqrt(secPtot2);                     // total P [GeV]
    double secEtot = vecCore::math::Sqrt(secPtot2 + secMass * secMass); // total energy in [GeV]
    double secEkin = secEtot - secMass;                  // kinetic energy in [GeV]
    // Ekin of the i-th secondary is higher than the threshold
    if (secEkin > energyLimit) { // insert secondary into tracks_v
      GeantTrack &track = td->GetNewTrack();
      // set the new track properties
      track.SetEvent(tracks.fEventV[iintrack]);
      track.SetEvslot(tracks.fEvslotV[iintrack]);
      track.SetPDG(secPDG);    // PDG code of this particle
      track.SetGVcode(pid[i]); // GV index of this particle
      track.SetCharge(secPartPDG->Charge()); // charge of this particle
#ifndef USE_VECGEOM_NAVIGATOR
      track.SetCharge(track.Charge()/3.);
#endif
      track.SetStatus(kNew);                  // status of this particle
      track.SetMass(secMass);                 // mass of this particle
      track.SetPosition(tracks.fXposV[iintrack], tracks.fYposV[iintrack], tracks.fZposV[iintrack]);
      track.SetDirection(px / secPtot, py / secPtot, pz / secPtot);
      track.SetP(secPtot);                    // momentum of this particle
      track.SetE(secEtot);                    // total E of this particle
      track.SetTime(tracks.fTimeV[iintrack]); // global time for this particle
      track.SetSafety(tracks.fSafetyV[iintrack]);
      track.SetBoundary(tracks.fBoundaryV[iintrack]);
      track.SetPath(tracks.fPathV[iintrack]);
      track.SetNextPath(tracks.fPathV[iintrack]);

      // rotate at-rest secondary by a common random theta and random phi
      RotateNewTrack(randDirX, randDirY, randDirZ, track);

      td->fPropagator->AddTrack(track);
      tracks.AddTrack(track);
      td->ReleaseTrack(track);

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
// will be called recursively if necessary
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::GetRestFinStates(int partindex, TMXsec *mxs, double energyLimit, GeantTrack *track,
                                   int &nTotSecPart, TrackVec_t &output, GeantTaskData *td) {
  // current track should have already been killed before calling
  const double mecc = 0.00051099906; // e- mass c2 in [GeV]
  double rndArray[3];
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
#ifdef USE_VECGEOM_NAVIGATOR
  td->fRndm->uniform_array(3, rndArray);
#elif USE_ROOT
  td->fRndm->RndmArray(3, rndArray);
#endif
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
    if (energyLimit > mecc) return;
    double randDirZ = 1.0 - 2.0 * rndArray[0];
    double randSinTheta = vecCore::math::Sqrt(1.0 - randDirZ * randDirZ);
    double randPhi = 2.0 * rndArray[1] * kPi;
    double randDirX = randSinTheta * vecCore::math::Cos(randPhi);
    double randDirY = randSinTheta * vecCore::math::Sin(randPhi);

    // need to do it one-by-one
    // 1. gamma
    GeantTrack &track1 = td->GetNewTrack();
    // set the new track properties: 2 gamma with m_{e}*c*c
    track1.SetEvent(track->Event());
    track1.SetEvslot(track->EventSlot());
    track1.SetPDG(22);                                    // gamma PDG code
    track1.SetGVcode(TPartIndex::I()->GetSpecGVIndex(2)); // gamma GV index
    track1.SetGeneration(track->GetGeneration() + 1);
    track1.SetStatus(kNew);                  // status of this particle
    track1.SetStage(int(kPreStepStage));
    track1.SetPosition(track->X(), track->Y(), track->Z());
    track1.SetDirection(randDirX, randDirY, randDirZ);
    track1.SetP(mecc);                       // momentum of this particle
    track1.SetE(mecc);                       // total E of this particle
    track1.SetTime(track->Time()); // total time of this particle
    track1.SetSafety(track->GetSafety());
    track1.SetBoundary(track->Boundary());
    track1.SetPath(track->Path());
    track1.SetNextPath(track->Path());

    td->fPropagator->AddTrack(track1);
    output.push_back(&track1);

    // 2. gamma : everything is the same but the direction
    GeantTrack &track2 = td->GetNewTrack();
    track2 = track1;
    track2.SetDirection(-randDirX, -randDirY, -randDirZ);

    td->fPropagator->AddTrack(track2);
    output.push_back(&track2);

    nTotSecPart += 2;
    return;
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

  // isSurv should always be false here because primary was stopped -> just a check
  if (isSurv)
    printf("A stopped particle survived its rest process in TTabPhysMgr::GetRestFinSates!\n");

  // for a random rotation
  if (nSecPart) {
    randDirZ = 1.0 - 2.0 * rndArray[1];
    randSinTheta = vecCore::math::Sqrt((1.0 - randDirZ) * (1.0 + randDirZ));
    randPhi = kTwoPi * rndArray[2];
    randDirX = randSinTheta * vecCore::math::Cos(randPhi);
    randDirY = randSinTheta * vecCore::math::Sin(randPhi);
  }

  // if tehere was any energy deposit add it to parent track doposited energy
  track->IncreaseEdep(kerma);

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
      track->IncreaseEdep(vecCore::math::Sqrt(secPtot2 + secMass * secMass) - secMass);
      continue;
    }

    int secPDG = TPartIndex::I()->PDG(pid[i]); // Geant V particle code -> particle PGD code
#ifdef USE_VECGEOM_NAVIGATOR
    const Particle_t *const &secPartPDG = &Particle_t::GetParticle(secPDG);
#else
    TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
#endif
    double secMass = secPartPDG->Mass();
    double px = mom[3 * i];
    double py = mom[3 * i + 1];
    double pz = mom[3 * i + 2];
    double secPtot2 = px * px + py * py + pz * pz;       // total P^2 [GeV^2]
    double secPtot = vecCore::math::Sqrt(secPtot2);                     // total P [GeV]
    double secEtot = vecCore::math::Sqrt(secPtot2 + secMass * secMass); // total energy in [GeV]
    double secEkin = secEtot - secMass;                  // kinetic energy in [GeV]
    // Ekin of the i-th secondary is higher than the threshold
    if (secEkin > energyLimit) { // insert secondary into tracks_v
      GeantTrack &track3 = td->GetNewTrack();
      // set the new track properties
      track3.SetEvent(track->Event());
      track3.SetEvslot(track->EventSlot());
      //       track3.fParticle = nTotSecPart;          //index of this particle
      track3.SetPDG(secPDG);    // PDG code of this particle
      track3.SetGVcode(pid[i]); // GV index of this particle
      track3.SetCharge(secPartPDG->Charge()); // charge of this particle
#ifndef USE_VECGEOM_NAVIGATOR
      track3.SetCharge(track3.Charge()/3.);
#endif
      track3.SetGeneration(track->GetGeneration() + 1);
      track3.SetStatus(kNew);                  // status of this particle
      track3.SetStage(int(kPreStepStage));
      track3.SetMass(secMass);                 // mass of this particle
      track3.SetPosition(track->X(), track->Y(), track->Z());
      track3.SetDirection(px / secPtot, py / secPtot, pz / secPtot);
      track3.SetP(secPtot);                    // momentum of this particle
      track3.SetE(secEtot);                    // total E of this particle
      track3.SetTime(track->Time()); // global time for this particle
      track3.SetSafety(track->GetSafety());
      track3.SetBoundary(track->Boundary());
      track3.SetPath(track->Path());
      track3.SetNextPath(track->Path());

      // rotate at-rest secondary by a common random theta and random phi
      RotateNewTrack(randDirX, randDirY, randDirZ, track3);

      td->fPropagator->AddTrack(track3);
      output.push_back(&track3);

      ++nTotSecPart; // increase # of secondaries in tracks_v
    } else {
      // add the Ekin of this secondary to the energy depositon
      track->IncreaseEdep(secEkin);
      // check if it is a stopped particle and call at-rest sampling if necessary
      if (secEkin <= 0.0 && HasRestProcess(pid[i]))
        GetRestFinStates(pid[i], mxs, energyLimit, track, nTotSecPart, output, td); // RECURSION
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
    Geant::Warning("TTabPhysMgr::SampleDecayInFlight","%s",
                   "\n---       A particle survived its decay!!!       ---\n"
                   "----    In TTabPhysMgr::SampleFinalStateAtRest     ---\n");

  if (nSecPart) {
    // Go for the secondaries
    double beta = tracks.fPV[iintrack] / tracks.fEV[iintrack];
    double bx = tracks.fXdirV[iintrack] * beta;
    double by = tracks.fYdirV[iintrack] * beta;
    double bz = tracks.fZdirV[iintrack] * beta;
    double b2 = bx * bx + by * by + bz * bz; // it is beta*beta
    double gam = 1.0 / vecCore::math::Sqrt(1.0 - b2);
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
        tracks.fEdepV[iintrack] += vecCore::math::Sqrt(secPtot2 + secMass * secMass) - secMass;
        continue;
      }

      int secPDG = TPartIndex::I()->PDG(pid[isec]); // GV part.code -> PGD code
#ifdef USE_VECGEOM_NAVIGATOR
      const Particle_t *const &secPartPDG = &Particle_t::GetParticle(secPDG);
#else
      TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
#endif
      double secMass = secPartPDG->Mass(); // mass [GeV]
      double px = mom[3 * isec];
      double py = mom[3 * isec + 1];
      double pz = mom[3 * isec + 2];
      double secP2 = px * px + py * py + pz * pz;       // total P^2 [GeV^2]
      double secEtot = vecCore::math::Sqrt(secP2 + secMass * secMass); // total E [GeV]
      // double secEkin  = secEtot - secMass; //kinetic energy in [GeV]

      double bp = bx * px + by * py + bz * pz;
      px = px + gam2 * bp * bx + gam * bx * secEtot;
      py = py + gam2 * bp * by + gam * by * secEtot;
      pz = pz + gam2 * bp * bz + gam * bz * secEtot;
      secEtot = gam * (secEtot + bp);

      double secPtot = vecCore::math::Sqrt((secEtot - secMass) * (secEtot + secMass));
      double secEkin = secEtot - secMass;
      if (secEkin > energyLimit) { // insert secondary into tracks_v
        GeantTrack &track = td->GetNewTrack();
        // set the new track properties
        track.SetEvent(tracks.fEventV[iintrack]);
        track.SetEvslot(tracks.fEvslotV[iintrack]);
        track.SetPDG(secPDG);       // PDG code of this particle
        track.SetGVcode(pid[isec]); // GV index of this particle
        track.SetCharge(secPartPDG->Charge()); // charge of this particle
#ifndef USE_VECGEOM_NAVIGATOR
	      track.SetCharge(track.Charge()/3.);
#endif
        track.SetStatus(kNew);                  // status of this particle
        track.SetMass(secMass);                 // mass of this particle
        track.SetPosition(tracks.fXposV[iintrack], tracks.fYposV[iintrack], tracks.fZposV[iintrack]);
        track.SetDirection(px / secPtot, py / secPtot, pz / secPtot);
        track.SetP(secPtot);                    // momentum of this particle
        track.SetE(secEtot);                    // total E of this particle
        track.SetTime(tracks.fTimeV[iintrack]); // global time for this track
        track.SetSafety(tracks.fSafetyV[iintrack]);
        track.SetBoundary(tracks.fBoundaryV[iintrack]);
        track.SetPath(tracks.fPathV[iintrack]);
        track.SetNextPath(tracks.fPathV[iintrack]);

        td->fPropagator->AddTrack(track);
        tracks.AddTrack(track);
        td->ReleaseTrack(track);

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

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::SampleDecayInFlight(int partindex, TMXsec *mxs, double energyLimit, GeantTrack *track,
                                      int &nTotSecPart, TrackVec_t &output, GeantTaskData *td) {
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
    double beta = track->Beta();
    double bx = track->Dx() * beta;
    double by = track->Dy() * beta;
    double bz = track->Dz() * beta;
    double b2 = bx * bx + by * by + bz * bz; // it is beta*beta
    double gam = 1.0 / vecCore::math::Sqrt(1.0 - b2);
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
        track->IncreaseEdep(vecCore::math::Sqrt(secPtot2 + secMass * secMass) - secMass);
        continue;
      }

      int secPDG = TPartIndex::I()->PDG(pid[isec]); // GV part.code -> PGD code
#ifdef USE_VECGEOM_NAVIGATOR
      const Particle_t *const &secPartPDG = &Particle_t::GetParticle(secPDG);
#else
      TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
#endif
      double secMass = secPartPDG->Mass(); // mass [GeV]
      double px = mom[3 * isec];
      double py = mom[3 * isec + 1];
      double pz = mom[3 * isec + 2];
      double secP2 = px * px + py * py + pz * pz;       // total P^2 [GeV^2]
      double secEtot = vecCore::math::Sqrt(secP2 + secMass * secMass); // total E [GeV]
      // double secEkin  = secEtot - secMass; //kinetic energy in [GeV]

      double bp = bx * px + by * py + bz * pz;
      px = px + gam2 * bp * bx + gam * bx * secEtot;
      py = py + gam2 * bp * by + gam * by * secEtot;
      pz = pz + gam2 * bp * bz + gam * bz * secEtot;
      secEtot = gam * (secEtot + bp);

      double secPtot = vecCore::math::Sqrt((secEtot - secMass) * (secEtot + secMass));
      double secEkin = secEtot - secMass;
      if (secEkin > energyLimit) { // insert secondary into tracks_v
        GeantTrack &track1 = td->GetNewTrack();
        // set the new track properties
        track1.SetEvent(track->Event());
        track1.SetEvslot(track->EventSlot());
        track1.SetPDG(secPDG);       // PDG code of this particle
        track1.SetGVcode(pid[isec]); // GV index of this particle
        track1.SetCharge(secPartPDG->Charge()); // charge of this particle
#ifndef USE_VECGEOM_NAVIGATOR
	      track1.SetCharge(track1.Charge()/3.);
#endif
        track1.SetGeneration(track->GetGeneration() + 1);
        track1.SetStatus(kNew);                  // status of this particle
        track1.SetStage(int(kPreStepStage));
        track1.SetMass(secMass);                 // mass of this particle
        track1.SetPosition(track->X(), track->Y(), track->Z());
        track1.SetDirection(px / secPtot, py / secPtot, pz / secPtot);
        track1.SetP(secPtot);                    // momentum of this particle
        track1.SetE(secEtot);                    // total E of this particle
        track1.SetTime(track->Time()); // global time for this track
        track1.SetSafety(track->GetSafety());
        track1.SetBoundary(track->Boundary());
        track1.SetPath(track->Path());
        track1.SetNextPath(track->Path());

        td->fPropagator->AddTrack(track1);
        output.push_back(&track1);

        ++nTotSecPart; // increase # of secondaries in tracks_v
      } else {
        // add the Ekin of this secondary to the energy depositon
        track->IncreaseEdep(secEkin);
        // check if it is a stopped particle and call at-rest sampling if necessary
        if (secEkin <= 0.0 && HasRestProcess(pid[isec]))
          GetRestFinStates(pid[isec], mxs, energyLimit, track, nTotSecPart, output, td);
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
  double sinTheta0 = vecCore::math::Sqrt(oldXdir * oldXdir + oldYdir * oldYdir);
  double cosPhi0;
  double sinPhi0;

  if (sinTheta0 > amin) {
    cosPhi0 = oldXdir / sinTheta0;
    sinPhi0 = oldYdir / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = track.Dx();
  double h1 = sinTheta0 * track.Dz() + cosTheta0 * h0;
  double h2 = track.Dy();

  track.SetDirection(h1 * cosPhi0 - h2 * sinPhi0,
                     h1 * sinPhi0 + h2 * cosPhi0,
                     track.Dz() * cosTheta0 - h0 * sinTheta0);

  // renormalization: -use 1-th order Taylor aprx. of 1/sqrt(x) around 1.0
  // that should be almost exact since the vector almost normalized!
  double delta = one5 - half * (track.Dx() * track.Dx() + track.Dy() * track.Dy() + track.Dz() * track.Dz());
  track.SetDirection(delta*track.Dx(), delta*track.Dy(), delta*track.Dz());
}

//______________________________________________________________________________
// FOR THE itrack-th element of a GeantTrack_v
// (oldXdir, oldYdir, oldZdir) is the direction vector of parent track in lab.
// frame; direction vector of the current track, measured from local Z is
// already updated in GeantTrack track; here we rotate it to lab. frame
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack_v &tracks, int itrack) {
  const double one = 1.0;
  const double zero = 0.0;
  const double amin = 1.0e-10;
  const double one5 = 1.5;
  const double half = 0.5;

  double cosTheta0 = oldZdir;
  double sinTheta0 = vecCore::math::Sqrt(oldXdir * oldXdir + oldYdir * oldYdir);
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
  double delta = one5 -
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

  double cosTheta0 = track.Dz();
  double sinTheta0 = vecCore::math::Sqrt(track.Dx() * track.Dx() + track.Dy() * track.Dy());
  double cosPhi0;
  double sinPhi0;
  double cosTheta = vecCore::math::Cos(theta);
  double sinTheta = vecCore::math::Sin(theta);

  if (sinTheta0 > amin) {
    cosPhi0 = track.Dx() / sinTheta0;
    sinPhi0 = track.Dy() / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = sinTheta * vecCore::math::Cos(phi);
  double h1 = sinTheta0 * cosTheta + cosTheta0 * h0;
  double h2 = sinTheta * vecCore::math::Sin(phi);

  track.SetDirection(h1 * cosPhi0 - h2 * sinPhi0,
                     h1 * sinPhi0 + h2 * cosPhi0,
                     cosTheta * cosTheta0 - h0 * sinTheta0);

  // renormalization: -ensure normality to avoid accumulated numerical errors
  //    due to sequential calls of rotation; avoid 1/sqrt(x) computation by
  //    using the 1-th order Taylor aprx. around 1.0 that should be almost
  //    exact since the vector almost normalized!
  double delta = one5 - half * (track.Dx() * track.Dx() + track.Dy() * track.Dy() + track.Dz() * track.Dz());
  track.SetDirection(delta * track.Dx(), delta * track.Dy(), delta * track.Dz());
}

//______________________________________________________________________________
// FOR THE itrack-th element of a GeantTrack_v
// GeantTrack_v contains the original direction in lab frame; theta and
// phi are the scattering angles measured form the particle local Z
VECCORE_ATT_HOST_DEVICE
void TTabPhysMgr::RotateTrack(GeantTrack_v &tracks, int itrack, double theta, double phi) {
  const double one = 1.0;
  const double zero = 0.0;
  const double amin = 1.0e-10;
  const double one5 = 1.5;
  const double half = 0.5;

  double cosTheta0 = tracks.fZdirV[itrack];
  double sinTheta0 =
      vecCore::math::Sqrt(tracks.fXdirV[itrack] * tracks.fXdirV[itrack] + tracks.fYdirV[itrack] * tracks.fYdirV[itrack]);
  double cosPhi0;
  double sinPhi0;
  double cosTheta = vecCore::math::Cos(theta);
  double sinTheta = vecCore::math::Sin(theta);

  if (sinTheta0 > amin) {
    cosPhi0 = tracks.fXdirV[itrack] / sinTheta0;
    sinPhi0 = tracks.fYdirV[itrack] / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = sinTheta * vecCore::math::Cos(phi);
  double h1 = sinTheta0 * cosTheta + cosTheta0 * h0;
  double h2 = sinTheta * vecCore::math::Sin(phi);

  tracks.fXdirV[itrack] = h1 * cosPhi0 - h2 * sinPhi0;
  tracks.fYdirV[itrack] = h1 * sinPhi0 + h2 * cosPhi0;
  tracks.fZdirV[itrack] = cosTheta * cosTheta0 - h0 * sinTheta0;

  // renormalization: -ensure normality to avoid accumulated numerical errors
  //    due to sequential calls of rotation; avoid 1/sqrt(x) computation by
  //    using the 1-th order Taylor aprx. around 1.0 that should be almost
  //    exact since the vector almost normalized!
  double delta = one5 -
                 half * (tracks.fXdirV[itrack] * tracks.fXdirV[itrack] + tracks.fYdirV[itrack] * tracks.fYdirV[itrack] +
                         tracks.fZdirV[itrack] * tracks.fZdirV[itrack]);
  tracks.fXdirV[itrack] *= delta;
  tracks.fYdirV[itrack] *= delta;
  tracks.fZdirV[itrack] *= delta;
}

//______________________________________________________________________________

const char *TTabPhysMgr::GetVersion() const {
#ifndef VECCORE_CUDA
  static bool first = true;
  static std::mutex l;
  static char ver[512];
  l.lock();
  if (first) {
    first = false;
    sprintf(ver, "%d.%d.%d", VersionMajor(), VersionMinor(), VersionSub());
  }
  l.unlock();
  return ver;
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool TTabPhysMgr::HasRestProcess(int gvindex) {
  return fDecay->HasDecay(gvindex) || fHasNCaptureAtRest[gvindex] || (gvindex == TPartIndex::I()->GetSpecGVIndex(1));
}

} // GEANT_IMPL_NAMESPACE
} // Geant
