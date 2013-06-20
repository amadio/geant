#ifndef GEANT_EMGAMMA
#define GEANT_EMGAMMA
// Simple physics processes for Gammas
// - Compton Scattering

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TGeoVolume;
class GeantTrack;
class ComptonCrossSection; 
class TRandom; 

#include "PhysicsProcess.h"
#include "TGeoMaterial.h"
#include "TLorentzVector.h" 

#include <vector>

//
//  Class for Compton Process
//  -------------------------
//  Duties:  to calculate cross section, and provide output of interaction if called
//
//   Implementation 1:  asks Geant4 process to calculate cross section, and interact
//
//   Author: John Apostolakis,  29 Feb 2012
//______________________________________________________________________________
class GammaCompton : public PhysicsProcess
{
public:
  GammaCompton(const char *name=""); 
  virtual ~GammaCompton() {}
  void SetTcut( Double_t tcut ) { fTcut= tcut; } 
  void SetTmax( Double_t tmax ) { fTmax= tmax; } 
  
  virtual void ComputeIntLen(TGeoVolume *vol,
                             Int_t ntracks, 
                             Int_t *trackin, 
                             Double_t *lengths, 
                             Int_t tid);
  virtual void PostStep(     TGeoVolume *vol,
                             Int_t ntracks, 
                             Int_t *trackin, 
                             Int_t &nout, 
                             Int_t* trackout, 
                             Int_t tid);

  void SampleSecondaries(    const TGeoMaterial&   tMaterial,
                             TLorentzVector&       gamma4Mom,     // In/Out: 4-mom of gamma
                             // Int_t&                 electronOut,   // Out: true if secondary created
                             TLorentzVector&       electron4Mom,  // Out: 4-mom of outgoing e-
                             // Double_t&             enDeposit,     // Out: Energy Deposit
                             TRandom*             fRndEngine); 
private:
  ComptonCrossSection   *fComptonXS; 
  Double_t               fTcut;
  Double_t               fTmax; 
  std::vector<Double_t>  fProbabilities; 

  ClassDef(GammaCompton,1)    // Compton process
};

#endif
