
#include "GSMSCModel.h"


// from material
#include "Types.h"
#include "Material.h"
#include "MaterialProperties.h"
#include "MaterialCuts.h"
#include "Element.h"

#include "Region.h"
#include "ELossTableManager.h"

#include "GSMSCTable.h"
#include "PWATotalXsecTable.h"

#include "Particle.h"
#include "Electron.h"
#include "Positron.h"

// from geantV
#include "GeantTaskData.h"
#include "GeantTrack.h"

#include <cmath>

namespace geantphysics {


GSMSCTable*        GSMSCModel::gGSTable      = nullptr;
PWATotalXsecTable* GSMSCModel::gPWAXsecTable = nullptr;


GSMSCModel::GSMSCModel(bool iselectron, const std::string &name) : fIsElectron(iselectron),  MSCModel(name) {
  fIsUseAccurate         = true;
  fIsOptimizationOn      = true;

  fIsUsePWATotalXsecData = false;

  fTauSmall       = 1.e-16;
  fTauLim         = 1.e-6;
  fTLimitMinfix2  = 1.*geant::nm;
  fDtrl           = 0.05;

  fParticle = Electron::Definition();
  if (!fIsElectron) {
    fParticle = Positron::Definition();
  }
  fCharge = fParticle->GetPDGCharge();
}


GSMSCModel::~GSMSCModel() {
  if (gGSTable) {
    delete gGSTable;
    gGSTable = nullptr;
  }
  if (gPWAXsecTable) {
    delete gPWAXsecTable;
    gPWAXsecTable = nullptr;
  }
}

void GSMSCModel::Initialize() {
  // create static GSTable and PWATotalXsecTable if needed
  if (!gGSTable) {
    gGSTable = new GSMSCTable();
  }
  // GSTable will be initialised (data are loaded) only if they are not there yet.
  gGSTable->Initialize();
  //
  if (fIsUsePWATotalXsecData) {
    if (!gPWAXsecTable) {
      gPWAXsecTable = new PWATotalXsecTable();
    }
    gPWAXsecTable->Initialise();
  }
}

void GSMSCModel::StepLimit(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) {
  bool   isOnBoundary         = gtrack->fBoundary;
  int    matIndx              = gtrack->GetMaterial()->GetIndex();
  int    regIndx              = const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetRegion()->GetIndex();
  const  MaterialCuts *matCut = MaterialCuts::GetMaterialCut(regIndx,matIndx);
  double kineticEnergy        = gtrack->fE-gtrack->fMass;

  if (kineticEnergy<GetLowEnergyUsageLimit() || kineticEnergy>GetHighEnergyUsageLimit()) {
    return;
  }

  double range                = ELossTableManager::Instance().GetRestrictedRange(matCut,fParticle,kineticEnergy);
  double skindepth            = 0.;
  double presafety            = gtrack->fSafety; // pre-step point safety
  double geomLimit            = gtrack->fSnext;  // init to distance-to-boundary
  //
  // Compute elastic mfp, first transport mfp, screening parameter and G1
  double lambel;
  double lambtr1;
  double scra;
  double g1;
  ComputeParameters(matCut, kineticEnergy, lambel, lambtr1, scra, g1);
  gtrack->fLambda0 = lambel;
  gtrack->fLambda1 = lambtr1;
  gtrack->fScrA    = scra;
  gtrack->fG1      = g1;
  gtrack->fRange   = range;
  //
  // Set initial values:
  //  : lengths are initialised to the current minimum physics step  which is the true, minimum
  //    step length from all other physics
  gtrack->fTheTrueStepLenght    = gtrack->fPstep;
  gtrack->fTheTransportDistance = gtrack->fPstep;
  gtrack->fTheZPathLenght       = gtrack->fPstep;
  gtrack->SetDisplacement(0.,0.,0.);
  gtrack->SetNewDirectionMsc(0.,0.,1.);

  // Can everything be done in the step limit phase ?
  gtrack->fIsEverythingWasDone  = false;
  // Multiple scattering needs to be sample ?
  gtrack->fIsMultipleSacettring = false;
  // Single scattering needs to be sample ?
  gtrack->fIsSingleScattering   = false;
  // Was zero deflection in multiple scattering sampling ?
  gtrack->fIsNoScatteringInMSC  = false;
  // Do not care about displacement in MSC sampling
  // ( used only in the case of fgIsOptimizationOn = true)
  gtrack->fIsNoDisplace         = false;

  // Zeff  = TotNbOfElectPerVolume/TotNbOfAtomsPerVolume
  double fZeff = matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()/
                 matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfAtomsPerVol();

  // distance will take into account max-fluct: we don't have fluctuation but use this for consistency
  double distance = range;
  distance *= (1.20-fZeff*(1.62e-2-9.22e-5*fZeff));
  //
  // Possible optimization : if the distance is samller than the safety -> the
  // particle will never leave this volume -> dispalcement
  // as the effect of multiple elastic scattering can be skipped
  // Important : this optimization can cause problems if one does scoring
  // in a bigger volume since MSC won't be done deep inside the volume when
  // distance < safety so don't use optimized-mode in such case.
  if (fIsOptimizationOn && (distance<presafety)) {
     // Indicate that we need to do MSC after transportation and no dispalcement.
     gtrack->fIsMultipleSacettring = true;
     gtrack->fIsNoDisplace = true;
  } else if (GetMSCSteppingAlgorithm()==MSCSteppingAlgorithm::kUseDistanceToBoundary) {
    //Compute geomlimit (and presafety) :
    // - geomlimit will be:
    //    == the straight line distance to the boundary if currentRange is
    //       longer than that
    //    == a big value [geombig = 1.e50*mm] if currentRange is shorter than
    //       the straight line distance to the boundary
    // - presafety will be updated as well
    // So the particle can travell 'gemlimit' distance (along a straight
    // line!) in its current direction:
    //  (1) before reaching a boundary (geomlimit < geombig) OR
    //  (2) before reaching its current range (geomlimit == geombig)
//      geomlimit = ComputeGeomLimit(track, presafety, currentRange);
    // Record that the particle is on a boundary
    if (isOnBoundary) {
      gtrack->fIsWasOnBoundary = true;
    }
//       if( (stepStatus==fGeomBoundary) || (stepStatus==fUndefined && presafety==0.0))
//          fIsWasOnBoundary = true;

    // Set skin depth = skin x elastic_mean_free_path
    skindepth = GetSkin()*lambel;
    // Init the flag that indicates that the particle is within a skindepth distance from a boundary
    gtrack->fIsInsideSkin = false;
    // Check if we can try Single Scattering because we are within skindepth
    // distance from/to a boundary OR the current minimum true-step-length is
    // shorter than skindepth. NOTICE: the latest has only efficieny reasons
    // because the MSC angular sampling is fine for any short steps but much
    // faster to try single scattering in case of short steps.
    if (isOnBoundary || (presafety<skindepth) || (gtrack->fTheTrueStepLenght<skindepth)) {
      // check if we are within skindepth distance from a boundary
      if (isOnBoundary || (presafety<skindepth)) {
        gtrack->fIsInsideSkin = true;
        gtrack->fIsWasOnBoundary = true; // NOTE:
      }
      //Try single scattering:
      // - sample distance to next single scattering interaction (sslimit)
      // - compare to current minimum length
      //      == if sslimit is the shorter:
      //          - set the step length to sslimit
      //          - indicate that single scattering needs to be done
      //      == else : nothing to do
      //- in both cases, the step length was very short so geometrical and
      //  true path length are the same
      double sslimit = -1.*lambel*std::log(td->fRndm->uniform());
      // compare to current minimum step length
      if (sslimit<gtrack->fTheTrueStepLenght) {
        gtrack->fTheTrueStepLenght  = sslimit;
        gtrack->fIsSingleScattering = true;
      }
      // short step -> true step length equal to geometrical path length
      gtrack->fTheZPathLenght  = gtrack->fTheTrueStepLenght;
      // Set taht everything is done in step-limit phase so no MSC call
      // We will check if we need to perform the single-scattering angular
      // sampling i.e. if single elastic scattering was the winer!
      gtrack->fIsEverythingWasDone = true;
      gtrack->fIsNoDisplace = true;  // NOTE:
    } else {
      // After checking we know that we cannot try single scattering so we will need to make an MSC step
      // Indicate that we need to make and MSC step. We do not check if we can do it now i.e. if
      // presafety>final_true_step_length so we let the fIsEverythingWasDone = false which indicates that
      // we will perform MSC after transportation.
      gtrack->fIsMultipleSacettring = true;
      // Init the first-real-step falg: it will indicate if we do the first
      // non-single scattering step in this volume with this particle
      gtrack->fIsFirstRealStep = false;
      // If previously the partcile was on boundary it was within skin as
      // well. When it is not within skin anymore it has just left the skin
      // so we make the first real MSC step with the particle.
      if (gtrack->fIsWasOnBoundary && !gtrack->fIsInsideSkin) {
        // reset the 'was on boundary' indicator flag
        gtrack->fIsWasOnBoundary = false;
        gtrack->fIsFirstRealStep = true;
      }
      // If this is the first-real msc step (the partcile has just left the skin) or this is the first step with
      // the particle (the particle has just born or primary):
      //   - set the initial range that will be used later to limit its step
      //     (only in this volume, because after boundary crossing at the
      //     first-real MSC step we will reset)
      //  - don't let the partcile to cross the volume just in one step
      if (gtrack->fIsFirstStep || gtrack->fIsFirstRealStep || gtrack->fTheInitialRange>1.e+20) { //NOTE:
        gtrack->fTheInitialRange = range;
        // If GeantTrack::fSnext(distance-to-boundary) < range then the particle might reach the boundary along its
        // initial direction before losing its energy (in this step). Otherwise, we can be sure that the particle will
        // lose its energy before reaching the boundary along a starigth line so there is no geometrical limit appalied.
        // [However, tgeom is set only in the first or the first-real MSC step. After the first or first real MSC step
        // the direction will change tgeom won't guaranty anything! But we will try to end up within skindepth from the
        // boundary using the actual value of distance-to-boundary(See later at step reduction close to boundary).]
        // fTrueGeomLimit this is the geometrical limit (distance-to-boundary) converted to true step length.
        if (geomLimit<range) {
          // transfrom straight line distance to the boundary to real step length based on the mean values (using the
          // prestep point first-transport mean free path i.e. no energy loss correction)
          if ((1.-geomLimit/lambtr1)>0.) {
            geomLimit = -lambtr1*std::log(1.-geomLimit/lambtr1);
          }
          // the 2-different cases that could lead us here:
          // 1.: gtrack->fIsFirstRealStep
          gtrack->fTheTrueGeomLimit = geomLimit;
          // 2.: gtrack->fIsFirstStep otherwise
          if (gtrack->fIsFirstStep) {
            gtrack->fTheTrueGeomLimit *= 2.;
          }
          gtrack->fTheTrueGeomLimit /= GetGeomFactor(); // facgeom
        } else {
          gtrack->fTheTrueGeomLimit = 1.e+20;
        }
      }
      // True step length limit from range factor. Noteice, that the initial range is used that was set at the first
      // step or first-real MSC step in this volume with this particle.
      double tlimit = GetRangeFactor()*gtrack->fTheInitialRange;
      // Take the minimum of the true step length limits coming from geometrical constraint or range-factor limitation
      tlimit = std::min(tlimit,gtrack->fTheTrueGeomLimit);
      // Step reduction close to boundary: we try to end up within skindepth from the boundary ( Notice: in case of
      // magnetic field it might not work because geomlimit is the straigth line distance to the boundary in the
      // currect direction and mag. field can change the initial direction. So the particle might hit some boundary
      // before in a different direction. However, here we restrict the true path length to this (straight line)
      // length so the corresponding transport distance (straight line) will be even shorter than
      // geomlimit-0.999*skindepth after the change of true->geom.
      if (geomLimit<range) {  // only if the particle can reach the boundary before losing its energy
         tlimit = std::min(tlimit, geomLimit-0.999*skindepth);
      }
      //
      // randomize 1st step or 1st 'normal' step in volume: according to Laszlo's algorithm
      if (gtrack->fIsFirstStep || gtrack->fIsFirstRealStep) {
        gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, RandomizeTrueStepLength(td, tlimit));
      } else {
        gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, tlimit);
      }
/*
      // randomize the true step length if msc limited the step
      if (tlimit<gtrack->fTheTrueStepLenght) {
        gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, RandomizeTrueStepLength(td, tlimit));
      } else {
      // just restrict the true step length to the msc limit
        gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, tlimit);
      }
*/
// NOTE: we changed this according to Urban ! See above the original one. (we probably should combine the 2 and use
//       this randomization only if msc limited the step in case of (1st or 1st normal steps))
/*
      // randomize 1st step or 1st 'normal' step in volume: according to Laszlo's algorithm
      if (gtrack->fIsFirstStep || gtrack->fIsFirstRealStep) {
        gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, RandomizeTrueStepLength());
      } else {
        gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, tlimit);
      }
*/
    }
  } else if (GetMSCSteppingAlgorithm()==MSCSteppingAlgorithm::kErrorFree) {
    geomLimit = presafety;  // pre-step point safety
    // Set skin depth = skin x elastic_mean_free_path
    skindepth = GetSkin()*lambel;
    // Check if we can try Single Scattering because we are within skindepth distance from/to a boundary OR the current
    // minimum true-step-length is shorter than skindepth. NOTICE: the latest has only efficieny reasons because the
    // MSC angular sampling is fine for any short steps but much faster to try single scattering in case of short steps.
    if (isOnBoundary || (presafety<skindepth) || (gtrack->fTheTrueStepLenght<skindepth)) {
      //Try single scattering:
      // - sample distance to next single scattering interaction (sslimit)
      // - compare to current minimum length
      //      == if sslimit is the shorter:
      //          - set the step length to sslimit
      //          - indicate that single scattering needs to be done
      //      == else : nothing to do
      //- in both cases, the step length was very short so geometrical and
      //  true path length are the same
      double sslimit = -1.*lambel*std::log(td->fRndm->uniform());
      // compare to current minimum step length
      if (sslimit<gtrack->fTheTrueStepLenght) {
        gtrack->fTheTrueStepLenght  = sslimit;
        gtrack->fIsSingleScattering = true;
      }
      // short step -> true step length equal to geometrical path length
      gtrack->fTheZPathLenght  = gtrack->fTheTrueStepLenght;
      // Set that everything is done in step-limit phase so no MSC call
      // We will check if we need to perform the single-scattering angular sampling i.e. if single elastic scattering
      // was the winer!
      gtrack->fIsEverythingWasDone = true;
      gtrack->fIsNoDisplace        = true; //NOTE:
    } else {
      // After checking we know that we cannot try single scattering so we will
      // need to make an MSC step
      // Indicate that we need to make and MSC step.
      gtrack->fIsMultipleSacettring = true;
      gtrack->fIsEverythingWasDone  = true;
      // limit from range factor
      gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, GetRangeFactor()*gtrack->fRange);
      // never let the particle go further than the safety if we are out of the skin
      // if we are here we are out of the skin, presafety > 0.
      if (gtrack->fTheTrueStepLenght>presafety) {
        gtrack->fTheTrueStepLenght = presafety;
      }
      // make sure that we are still within the aplicability of condensed histry model
      // i.e. true step length is not longer than first transport mean free path.
      // We schould take into account energy loss along 0.5x lambda_transport1
      // step length as well. So let it 0.4 x lambda_transport1
      gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, 0.4*lambtr1);
    }
  } else if (GetMSCSteppingAlgorithm()==MSCSteppingAlgorithm::kUseSaftey) {
    // This is the default stepping algorithm: the fastest but the least
    // accurate that corresponds to fUseSafety in Urban model. Note, that GS
    // model can handle any short steps so we do not need the minimum limits

    // NO single scattering in case of skin or short steps (by defult the MSC
    // model will be single or even no scattering in case of short steps
    // compared to the elastic mean free path.)

    // indicate that MSC needs to be done (always and always after transportation)
    gtrack->fIsMultipleSacettring = true;
    // far from boundary-> in optimized mode do not sample dispalcement. (safety is zero if we are on boundary)
    if ((distance < presafety) && (fIsOptimizationOn)) {
      gtrack->fIsNoDisplace = true;
    } else {
      // Urban like
      double fr = GetRangeFactor();
      if (gtrack->fIsFirstStep || isOnBoundary || gtrack->fTheInitialRange>1.e+20) { //NOTE:
        gtrack->fTheInitialRange = range;
// We don't use this: we won't converge to the single scattering results with
//                    decreasing range-factor.
//              rangeinit = std::max(rangeinit, fLambda1);
//              if(fLambda1 > lambdalimit) {
//                fr *= (0.75+0.25*fLambda1/lambdalimit);
//              }
      }
      //step limit
      double tlimit = std::max(fr*gtrack->fTheInitialRange, GetSafetyFactor()*presafety);
      // first step randomization
      if (gtrack->fIsFirstStep || isOnBoundary) { // NOTE: we should probably randomize only if the step was limited by msc
        gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, RandomizeTrueStepLength(td, tlimit));
      } else {
        gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght, tlimit);
      }
    }
  }
  // msc step limit is done!
  //
  //
  // reset first step flag
  gtrack->fIsFirstStep =false;
  // performe single scattering, multiple scattering if this later can be done safely here
  if (gtrack->fIsEverythingWasDone) {
    if (gtrack->fIsSingleScattering) {
      // sample single scattering
      double cost,sint,cosPhi,sinPhi;
      SingleScattering(scra, td, cost, sint);
      double phi = geant::kTwoPi*td->fRndm->uniform();
      sinPhi = std::sin(phi);
      cosPhi = std::cos(phi);
      gtrack->SetNewDirectionMsc(sint*cosPhi,sint*sinPhi,cost);
    } else if (gtrack->fIsMultipleSacettring) {
      // sample multiple scattering
      SampleMSC(gtrack, td); // fTheZPathLenght, fTheDisplacementVector and fTheNewDirection will be set
    } // and if single scattering but it was longer => nothing to do
  } //else { do nothing here but after transportation
  //
  // convert the true step length to geometrical one
  ConvertTrueToGeometricLength(gtrack, td);
}

bool GSMSCModel::SampleScattering(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) {
  if (GetMSCSteppingAlgorithm()==MSCSteppingAlgorithm::kUseDistanceToBoundary && gtrack->fIsEverythingWasDone && gtrack->fIsSingleScattering) { // ONLY single scattering is done in advance
    // single scattering was and scattering happend
    RotateToLabFrame(gtrack->fTheNewDirectionX, gtrack->fTheNewDirectionY, gtrack->fTheNewDirectionZ,
                     gtrack->fXdir, gtrack->fYdir, gtrack->fZdir);
    // displacement is left to (0,0,0)
    //fParticleChange->ProposeMomentumDirection(fTheNewDirection);
    return true; // i.e. new direction
  } else if (GetMSCSteppingAlgorithm()==MSCSteppingAlgorithm::kErrorFree) {
    if (gtrack->fIsEndedUpOnBoundary) {// do nothing
      // displacement is left to (0,0,0)
      return false; // i.e. no new direction
    } else if (gtrack->fIsEverythingWasDone) { // evrything is done if not optimizations case !!!
      // check single scattering and see if it happened
      if (gtrack->fIsSingleScattering) {
        RotateToLabFrame(gtrack->fTheNewDirectionX, gtrack->fTheNewDirectionY, gtrack->fTheNewDirectionZ,
                         gtrack->fXdir, gtrack->fYdir, gtrack->fZdir);
        // displacement is left to (0,0,0)
        //fParticleChange->ProposeMomentumDirection(fTheNewDirection);
        return true; // i.e. new direction
      }
      // check if multiple scattering happened and do things only if scattering was really happening
      if (gtrack->fIsMultipleSacettring && !gtrack->fIsNoScatteringInMSC) {
        RotateToLabFrame(gtrack->fTheNewDirectionX, gtrack->fTheNewDirectionY, gtrack->fTheNewDirectionZ,
                         gtrack->fXdir, gtrack->fYdir, gtrack->fZdir);
        RotateToLabFrame(gtrack->fTheDisplacementVectorX, gtrack->fTheDisplacementVectorY, gtrack->fTheDisplacementVectorZ,
                         gtrack->fXdir, gtrack->fYdir, gtrack->fZdir);
        return true; // i.e. new direction
      }
      // The only thing that could happen if we are here (kErrorFree and fIsEverythingWasDone)
      // is that  single scattering was tried but did not win so scattering did not happen.
      // So no displacement and no scattering
      return false; // i.e. no new direction
    }
    //
    // The only thing that could still happen with kErrorFree is that we are in the
    // optimization branch: so sample MSC angle here (no displacement)
  }
  // else MSC needs to be done here
  SampleMSC(gtrack, td);
  if (!gtrack->fIsNoScatteringInMSC) {
    RotateToLabFrame(gtrack->fTheNewDirectionX, gtrack->fTheNewDirectionY, gtrack->fTheNewDirectionZ,
                     gtrack->fXdir, gtrack->fYdir, gtrack->fZdir);
    if (!gtrack->fIsNoDisplace) {
      RotateToLabFrame(gtrack->fTheDisplacementVectorX, gtrack->fTheDisplacementVectorY, gtrack->fTheDisplacementVectorZ,
                       gtrack->fXdir, gtrack->fYdir, gtrack->fZdir);
    }
    return true; // new direction
  }
  return false;
}


void GSMSCModel::ConvertTrueToGeometricLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) {
  gtrack->fPar1 = -1.;
  gtrack->fPar2 =  0.;
  gtrack->fPar3 =  0.;
  // if fIsEverythingWasDone = TRUE  => fTheZPathLenght is already set so return with the already known value
  // Otherwise:
  if (!gtrack->fIsEverythingWasDone) {
    // this correction needed to run MSC with eIoni and eBrem inactivated
    // and makes no harm for a normal run
    gtrack->fTheTrueStepLenght = std::min(gtrack->fTheTrueStepLenght,gtrack->fRange);
    //  do the true -> geom transformation
    gtrack->fTheZPathLenght = gtrack->fTheTrueStepLenght;
    // z = t for very small true-path-length
    if (gtrack->fTheTrueStepLenght<fTLimitMinfix2) {
      return;
    }
    //
    double ekin = gtrack->fE-gtrack->fMass;
    double tau  = gtrack->fTheTrueStepLenght/gtrack->fLambda1;
    if (tau<=fTauSmall) {
      gtrack->fTheZPathLenght = std::min(gtrack->fTheTrueStepLenght, gtrack->fLambda1);
    } else if (gtrack->fTheTrueStepLenght<gtrack->fRange*fDtrl) {
      if (tau<fTauLim) {
        gtrack->fTheZPathLenght = gtrack->fTheTrueStepLenght*(1.-0.5*tau) ;
      } else {
        gtrack->fTheZPathLenght = gtrack->fLambda1*(1.-std::exp(-tau));
      }
    } else if(ekin<geant::kElectronMassC2 || gtrack->fTheTrueStepLenght==gtrack->fRange)  {
      gtrack->fPar1 = 1./gtrack->fRange ;                    // alpha =1/range_init for Ekin<mass
      gtrack->fPar2 = 1./(gtrack->fPar1*gtrack->fLambda1) ;  // 1/(alphaxlambda01)
      gtrack->fPar3 = 1.+gtrack->fPar2 ;
      if (gtrack->fTheTrueStepLenght<gtrack->fRange) {
        gtrack->fTheZPathLenght = 1./(gtrack->fPar1*gtrack->fPar3)
                                  * (1.-std::pow(1.-gtrack->fPar1*gtrack->fTheTrueStepLenght,gtrack->fPar3));
      } else {
        gtrack->fTheZPathLenght = 1./(gtrack->fPar1*gtrack->fPar3);
      }
    } else {
      int    matIndx              = gtrack->GetMaterial()->GetIndex();
      int    regIndx              = const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetRegion()->GetIndex();
      const  MaterialCuts *matCut = MaterialCuts::GetMaterialCut(regIndx,matIndx);
      double rfin    = std::max(gtrack->fRange-gtrack->fTheTrueStepLenght, 0.01*gtrack->fRange);
      double T1      = ELossTableManager::Instance().GetEnergyForRestrictedRange(matCut,fParticle,rfin);
      double lambda1 = GetTransportMeanFreePathOnly(matCut,T1);
      gtrack->fPar1  = (gtrack->fLambda1-lambda1)/(gtrack->fLambda1*gtrack->fTheTrueStepLenght);  // alpha
      gtrack->fPar2  = 1./(gtrack->fPar1*gtrack->fLambda1);
      gtrack->fPar3  = 1.+gtrack->fPar2;
      gtrack->fTheZPathLenght = 1./(gtrack->fPar1*gtrack->fPar3)
                                * (1.-std::pow(1.-gtrack->fPar1*gtrack->fTheTrueStepLenght,gtrack->fPar3));
    }
  }
  gtrack->fTheZPathLenght = std::min(gtrack->fTheZPathLenght,gtrack->fLambda1); // NOTE:
}


void GSMSCModel::ConvertGeometricToTrueLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) {
  // init
  gtrack->fIsEndedUpOnBoundary = false;
  // step was not defined by transportation: i.e. physics
  if (!gtrack->fBoundary) {
//  if ( std::abs(gtrack->fStep-gtrack->fTheZPathLenght)<1.e-8) {
    return; // fTheTrueStepLenght is known because the particle went as far as we expected
  }
  // else ::
  // - set the flag that transportation was the winer so DoNothin in DOIT !!
  // - convert geom -> true by using the mean value
  gtrack->fIsEndedUpOnBoundary = true; // OR LAST STEP
  // get the geometrical step length
  gtrack->fTheZPathLenght      = gtrack->fStep;
  // was a short single scattering step
  if (gtrack->fIsEverythingWasDone && !gtrack->fIsMultipleSacettring) {
    gtrack->fTheTrueStepLenght = gtrack->fStep;
    return;
  }
  // t = z for very small step
  if (gtrack->fStep<fTLimitMinfix2) {
    gtrack->fTheTrueStepLenght = gtrack->fStep;
    // recalculation
  } else {
    double tlength = gtrack->fStep;
    if (gtrack->fStep>gtrack->fLambda1*fTauSmall) {
      if (gtrack->fPar1< 0.) {
        tlength = -gtrack->fLambda1*std::log(1.-gtrack->fStep/gtrack->fLambda1) ;
      } else {
        double dum = gtrack->fPar1*gtrack->fPar3*gtrack->fStep;
        if (dum<1.) {
          tlength = (1.-std::pow(1.-dum,1./gtrack->fPar3))/gtrack->fPar1;
        } else {
          tlength = gtrack->fRange;
        }
      }
      if (tlength<gtrack->fStep || tlength>gtrack->fTheTrueStepLenght) {
        tlength = gtrack->fStep;
      }
    }
    gtrack->fTheTrueStepLenght = tlength;
  }
}

void GSMSCModel::SampleMSC(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) {
  gtrack->fIsNoScatteringInMSC = false;
  //
  int    matIndx              = gtrack->GetMaterial()->GetIndex();
  int    regIndx              = const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetRegion()->GetIndex();
  const  MaterialCuts *matCut = MaterialCuts::GetMaterialCut(regIndx,matIndx);
  double kineticEnergy        = gtrack->fE-gtrack->fMass;
  double range                = gtrack->fRange;             // set in the step limit phase
  double trueStepL            = gtrack->fTheTrueStepLenght; // proposed by all other physics
  //
  // Energy loss correction (2 versions): accurate and approximate
  //
  // what will be the energy loss due to along step energy losses if the particle goes to the current physics step
  double eloss  = kineticEnergy-ELossTableManager::Instance().GetEnergyForRestrictedRange(matCut,fParticle,range-trueStepL);

//  if (fTheTrueStepLenght > currentRange*dtrl) {
//    eloss = kineticEnergy-
//      GetEnergy(particle,range-gtrack->fTheTrueStepLenght,currentCouple);
//  } else {
//    eloss = fTheTrueStepLenght*GetDEDX(particle,kineticEnergy,currentCouple);
//  }

  double tau  = 0.;//    = kineticEnergy/electron_mass_c2; // where kinEnergy is the mean kinetic energy
  double tau2 = 0.;//   = tau*tau;
  double eps0 = 0.;//   = eloss/kineticEnergy0; // energy loss fraction to the begin step energy
  double epsm = 0.;//   = eloss/kineticEnergy;  // energy loss fraction to the mean step energy
  //
  // - init.
  double efEnergy = kineticEnergy;
  double efStep   = trueStepL;
  //
  double kineticEnergy0 = kineticEnergy;
  if (fIsUseAccurate) {    // - use accurate energy loss correction
    kineticEnergy -= 0.5*eloss;  // mean energy along the full step
    // other parameters for energy loss corrections
    tau    = kineticEnergy/geant::kElectronMassC2; // where kinEnergy is the mean kinetic energy
    tau2   = tau*tau;
    eps0   = eloss/kineticEnergy0; // energy loss fraction to the begin step energy
    epsm   = eloss/kineticEnergy;  // energy loss fraction to the mean step energy

    efEnergy    = kineticEnergy * (1. - epsm*epsm*(6.+10.*tau+5.*tau2)/(24.*tau2+48.*tau+72.));
    double dum  = 0.166666*(4.+tau*(6.+tau*(7.+tau*(4.+tau))))*(epsm/((tau+1.)*(tau+2.)))*(epsm/((tau+1.)*(tau+2.)));
    efStep = trueStepL*(1.-dum);
  } else {                  // - take only mean energy
      kineticEnergy  -= 0.5*eloss;  // mean energy along the full step
      efEnergy        = kineticEnergy;
      double factor   = 1./(1. + kineticEnergy/(2.*geant::kElectronMassC2));
      eps0            = eloss/kineticEnergy0;
      epsm            = eps0/(1.-0.5*eps0);
      double dum      = 0.3*(1.-factor*(1.-0.333333*factor))*eps0*eps0;
      efStep          = trueStepL*(1.+dum);
  }
  //
  // Get elastic mfp, first transport mfp, screening parameter, and G1 (were computed at pre-step step limit and stored)
  double lambel;
  double lambtr1;
  double scra;
  double g1;
  ComputeParameters(matCut, efEnergy, lambel, lambtr1, scra, g1);
  // s/lambda_el i.e. mean number of elastic scattering along the step
  double lambdan=0.;
  if (lambel>0.0) {
    lambdan = efStep/lambel;
  }
  if (lambdan<=1.0e-12) {  //
    if (gtrack->fIsEverythingWasDone) {
      gtrack->fTheZPathLenght = trueStepL;
    }
    gtrack->fIsNoScatteringInMSC = true;
    return;
  }
  //
  // correction from neglected term 1+A (only for Moliere parameter) NOTE: should correct lambtr1 as well because it
  // depends on lambel i.e. lambtr1=lambel/g1 !!!
  if (!fIsUsePWATotalXsecData) {
    lambdan = lambdan/(1.+scra);
  }
  double Qn1 = lambdan*g1; //2.* lambdan *scrA*((1.+scrA)*log(1.+1./scrA)-1.);
  //
  // Sample scattering angles
  // new direction, relative to the orriginal one is in {uss,vss,wss}
  double cosTheta1 = 1.0, sinTheta1 = 0.0, cosTheta2 = 1.0, sinTheta2 = 0.0;
  double cosPhi1   = 1.0, sinPhi1   = 0.0, cosPhi2   = 1.0, sinPhi2   = 0.0;
  double uss = 0.0, vss = 0.0, wss = 1.0;
  double x_coord = 0.0, y_coord = 0.0, z_coord = 1.0;
  double u2 = 0.0, v2 = 0.0;
  // if we are above the upper grid limit with lambdaxG1=true-length/first-trans-mfp
  // => izotropic distribution: lambG1_max =7.99 but set it to 7
  double *rndArray = td->fDblArray;
  if (0.5*Qn1>7.0) {
    //get 2 random numbers
    td->fRndm->uniform_array(2, rndArray);
    cosTheta1 = 1.-2.*rndArray[0];
    sinTheta1 = std::sqrt((1.-cosTheta1)*(1.+cosTheta1));
    cosTheta2 = 1.-2.*rndArray[1];
    sinTheta2 = std::sqrt((1.-cosTheta2)*(1.+cosTheta2));
  } else {
     // sample 2 scattering cost1, sint1, cost2 and sint2 for half path
     gGSTable->Sampling(0.5*lambdan, 0.5*Qn1, scra, cosTheta1, sinTheta1, td);
     gGSTable->Sampling(0.5*lambdan, 0.5*Qn1, scra, cosTheta2, sinTheta2, td);
     if (cosTheta1+cosTheta2>=2.) { // no scattering happened
        if (gtrack->fIsEverythingWasDone) {
           gtrack->fTheZPathLenght = trueStepL;
        }
        gtrack->fIsNoScatteringInMSC = true;
        return;
     }
  }
  // sample 2 azimuthal angles
  // get 2 random numbers
  td->fRndm->uniform_array(2, rndArray);
  double phi1 = geant::kTwoPi*rndArray[0];
  sinPhi1     = std::sin(phi1);
  cosPhi1     = std::cos(phi1);
  double phi2 = geant::kTwoPi*rndArray[1];
  sinPhi2     = std::sin(phi2);
  cosPhi2     = std::cos(phi2);
  // compute final direction realtive to z-dir
  u2  = sinTheta2*cosPhi2;
  v2  = sinTheta2*sinPhi2;
  double u2p = cosTheta1*u2 + sinTheta1*cosTheta2;
  uss  = u2p*cosPhi1 - v2*sinPhi1;
  vss  = u2p*sinPhi1 + v2*cosPhi1;
  wss  = cosTheta1*cosTheta2 - sinTheta1*u2;
  //
  // set the new direction proposed by msc: will be applied if the step doesn't end on boundary
  gtrack->SetNewDirectionMsc(uss,vss,wss);
  // set the fTheZPathLenght if we don't sample displacement and
  // we should do everything at the step-limit-phase before we return
  if (gtrack->fIsNoDisplace && gtrack->fIsEverythingWasDone) {
    gtrack->fTheZPathLenght = trueStepL;
  }
  //
  // in optimized-mode if the current-safety > current-range we do not use dispalcement
  if (gtrack->fIsNoDisplace) {
    return;
  }
  //
  //
  // Compute final position
  if (fIsUseAccurate) {
    // correction parameter
    double par = 1.;
    if (Qn1<0.7) {
      par = 1.;
    } else if (Qn1<7.0) {
      par = -0.031376*Qn1+1.01356;
    } else {
      par = 0.79;
    }
    // Moments with energy loss correction
    // --first the uncorrected (for energy loss) values of gamma, eta, a1=a2=0.5*(1-eta), delta
    // gamma = G_2/G_1 based on G2 computed from A by using the Wentzel DCS form of G2
    double loga   = std::log(1.0+1.0/scra);
    double gamma  = 6.0*scra*(1.0+scra)*(loga*(1.0+2.0*scra)-2.0)/g1;
    // sample eta from p(eta)=2*eta i.e. P(eta) = eta_square ;-> P(eta) = rand --> eta = sqrt(rand)
    double eta    = std::sqrt(td->fRndm->uniform());
    double eta1   = 0.5*(1-eta);  // used  more than once
    // 0.5 +sqrt(6)/6 = 0.9082483;
    // 1/(4*sqrt(6))  = 0.1020621;
    // (4-sqrt(6)/(24*sqrt(6))) = 0.026374715
    // delta = 0.9082483-(0.1020621-0.0263747*gamma)*Qn1 without energy loss cor.
    double delta  = 0.9082483-(0.1020621-0.0263747*gamma)*Qn1;
    //
    // compute alpha1 and alpha2 for energy loss correction
    double temp1  = 2.0+tau;
    double temp   = (2.0+tau*temp1)/((tau+1.0)*temp1);
    //take into account logarithmic dependence
    temp  = temp - (tau+1.0)/((tau+2.0)*(loga*(1.0+scra)-1.0));
    temp  = temp * epsm;
    temp1 = 1.0 - temp;
    delta = delta + 0.40824829*(eps0*(tau+1.0)/((tau+2.0)*(loga*(1.0+scra)-1.0)*(loga*(1.0+2.0*scra)-2.0)) - 0.25*temp*temp);
    double b = eta*delta;
    double c = eta*(1.0-delta);
    //
    // calculate transport direction cosines:
    // ut,vt,wt is the final position divided by the true step length
    double w1v2 = cosTheta1*v2;
    double ut   = b*sinTheta1*cosPhi1 + c*(cosPhi1*u2 - sinPhi1*w1v2) + eta1*uss*temp1;
    double vt   = b*sinTheta1*sinPhi1 + c*(sinPhi1*u2 + cosPhi1*w1v2) + eta1*vss*temp1;
    double wt   = eta1*(1+temp) +       b*cosTheta1 +  c*cosTheta2    + eta1*wss*temp1;
    //
    // long step correction (only if not error-free stepping algorithm is used)
    ut *=par;
    vt *=par;
    wt *=par;
    //
    // final position relative to the pre-step point in the scattering frame
    // ut = x_f/s so needs to multiply by s
    x_coord = ut*trueStepL;
    y_coord = vt*trueStepL;
    z_coord = wt*trueStepL;
    //
    if (gtrack->fIsEverythingWasDone) {
      // We sample in the step limit so set fTheZPathLenght = transportDistance
      // and lateral displacement (x_coord,y_coord,z_coord-transportDistance)
      // Calculate transport distance
      double transportDistance  = std::sqrt(x_coord*x_coord+y_coord*y_coord+z_coord*z_coord);
      // protection
      if (transportDistance>trueStepL) {
        transportDistance = trueStepL;
      }
      gtrack->fTheZPathLenght = transportDistance;
    }
    // else:: we sample in the post step point so
    //       the fTheZPathLenght was already set and was taken as transport along zet
    gtrack->SetDisplacement(x_coord,y_coord,z_coord-gtrack->fTheZPathLenght);
  } else {
    // compute zz = <z>/tPathLength
    // s -> true-path-length
    // z -> geom-path-length:: when PRESTA is used z =(def.) <z>
    // r -> lateral displacement = s/2 sin(theta)  => x_f = r cos(phi); y_f = r sin(phi)
    double zz = 0.0;
    if (gtrack->fIsEverythingWasDone) {
      // We sample in the step limit so set fTheZPathLenght = transportDistance
      // and lateral displacement (x_coord,y_coord,z_coord-transportDistance)
      if (Qn1<0.1) { // use 3-order Taylor approximation of (1-exp(-x))/x around x=0
        zz = 1.0 - Qn1*(0.5 - Qn1*(0.166666667 - 0.041666667*Qn1)); // 1/6 =0.166..7 ; 1/24=0.041..
      } else {
        zz = (1.-std::exp(-Qn1))/Qn1;
      }
    } else {
      // we sample in the post step point so
      // the fTheZPathLenght was already set and was taken as transport along zet
      zz = gtrack->fTheZPathLenght/trueStepL;
    }
    //
    double rr = (1.-zz*zz)/(1.-wss*wss); // s^2 >= <z>^2+r^2  :: where r^2 = s^2/4 sin^2(theta)
    if (rr>=0.25) {  // (1-<z>^2/s^2)/sin^2(theta) >= r^2/(s^2 sin^2(theta)) = 1/4 must hold
      rr = 0.25;
    }
    double rperp = trueStepL*std::sqrt(rr);  // this is r/sint
    x_coord  = rperp*uss;
    y_coord  = rperp*vss;
    z_coord  = zz*trueStepL;
    if (gtrack->fIsEverythingWasDone) {
      double transportDistance = std::sqrt(x_coord*x_coord + y_coord*y_coord + z_coord*z_coord);
      gtrack->fTheZPathLenght  = transportDistance;
    }
    gtrack->SetDisplacement(x_coord,y_coord,z_coord-gtrack->fTheZPathLenght);
  }
}


//
// sample single scattering based on the Wentzel DCS
void GSMSCModel::SingleScattering(double scra, Geant::GeantTaskData *td, double &cost, double &sint) {
  double rand1 = td->fRndm->uniform();
  // sampling 1-cos(theta)
  double dum0  = 2.0*scra*rand1/(1.0 -rand1+scra);
  // add protections
  if (dum0<0.0) {
    dum0 = 0.0;
  }  else if (dum0>2.0 ) {
    dum0 = 2.0;
  }
  // compute cos(theta) and sin(theta)
  cost = 1.0-dum0;
  sint = std::sqrt(dum0*(2.0-dum0));
}


//
// computes the elastic and first transport mfp, the screening parameter and the first transport coeficient values
//GetTransportMeanFreePath
void GSMSCModel::ComputeParameters(const MaterialCuts *matcut, double ekin, double &lambel, double &lambtr1,
                                   double &scra, double &g1) {
  double efEnergy = ekin;
  // NOTE: change this !!! changed !
  if (efEnergy<GetLowEnergyUsageLimit() ) {
    efEnergy = GetLowEnergyUsageLimit();
  }
  if (efEnergy>GetHighEnergyUsageLimit()) {
    efEnergy = GetHighEnergyUsageLimit();
  }
  //
  const Material*  theMaterial = matcut->GetMaterial();
  //
  lambel  = 0.0; // elastic mean free path
  lambtr1 = 0.0; // first transport mean free path
  scra    = 0.0; // screening parameter
  g1      = 0.0; // first transport coef.
  //
  // NOTE: change this !!! : add protection that do not use pwa screening above "a given energy ?": changed
  if (fIsUsePWATotalXsecData && efEnergy<PWATotalXsecZ::GetHighestEnergy()) {
    // use PWA total xsec data to determine screening and lambda0 lambda1
    const Vector_t<Element*> &theElemVect   = theMaterial->GetElementVector();
    const double* theAtomicNumDensityVector = theMaterial->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
    int elowindx = gPWAXsecTable->GetPWATotalXsecForZet(std::lrint(theElemVect[0]->GetZ()))
                   ->GetPWATotalXsecEnergyBinIndex(efEnergy);
    int numElems = theMaterial->GetNumberOfElements();
    for (int i=0; i<numElems; ++i) {
      int zet       = std::lrint(theElemVect[i]->GetZ());
      // total elastic scattering cross section in Geant4 internal length2 units
      int indx      = (int)(1.5+fCharge*1.5); // specify that want to get the total elastic scattering xsection
      double elxsec = (gPWAXsecTable->GetPWATotalXsecForZet(zet))->GetInterpXsec(efEnergy,elowindx,indx);
      // first transport cross section in Geant4 internal length2 units
      indx          = (int)(2.5+fCharge*1.5); // specify that want to get the first tarnsport xsection
      double t1xsec = (gPWAXsecTable->GetPWATotalXsecForZet(zet))->GetInterpXsec(efEnergy,elowindx,indx);
      lambel  += (theAtomicNumDensityVector[i]*elxsec);
      lambtr1 += (theAtomicNumDensityVector[i]*t1xsec);
    }
    if (lambel>0.0) {
      lambel =1./lambel;
    }
    if (lambtr1>0.0) {
      lambtr1 = 1./lambtr1;
      g1      = lambel/lambtr1;
    }
    scra = gGSTable->GetScreeningParam(g1);
  } else {
    // Get SCREENING FROM MOLIER
    // below 1 keV it can give bananas so prevet (1 keV)
    if (efEnergy<1.*geant::keV) {
       efEnergy = 1.*geant::keV;
    }
    // total mometum square in geantV internal energy2 units which is GeV2
    double pt2     = efEnergy*(efEnergy+2.0*geant::kElectronMassC2);
    double beta2   = pt2/(pt2+geant::kElectronMassC2*geant::kElectronMassC2);
    int    matindx = theMaterial->GetIndex();
    double bc      = gGSTable->GetMoliereBc(matindx);
    scra           = gGSTable->GetMoliereXc2(matindx)/(4.0*pt2*bc);
    // total elastic mean free path in geantV internal lenght units
    lambel         = beta2/bc;//*(1.+fScrA);  //NOTE: we should apply correction to the lambtr1 as well !
    g1             = 2.0*scra*((1.0+scra)*std::log(1.0/scra+1.0)-1.0);
    lambtr1        = lambel/g1;
  }
}


double GSMSCModel::GetTransportMeanFreePathOnly(const MaterialCuts *matcut, double ekin) {
  double efEnergy = ekin;
  // NOTE: change this !!! changed !
  if (efEnergy<GetLowEnergyUsageLimit() ) {
    efEnergy = GetLowEnergyUsageLimit();
  }
  if (efEnergy>GetHighEnergyUsageLimit()) {
    efEnergy = GetHighEnergyUsageLimit();
  }
  //
  const Material*  theMaterial = matcut->GetMaterial();
  //
  double lambtr1 = 0.0; // first transport mean free path
  if (fIsUsePWATotalXsecData && efEnergy<PWATotalXsecZ::GetHighestEnergy()) {
    // use PWA total xsec data to determine screening and lambda0 lambda1
    const Vector_t<Element*> &theElemVect   = theMaterial->GetElementVector();
    const double* theAtomicNumDensityVector = theMaterial->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
    int elowindx = gPWAXsecTable->GetPWATotalXsecForZet(std::lrint(theElemVect[0]->GetZ()))
                   ->GetPWATotalXsecEnergyBinIndex(efEnergy);
    int numElems = theMaterial->GetNumberOfElements();
    for (int i=0; i<numElems; ++i) {
      int zet       = std::lrint(theElemVect[i]->GetZ());
      // first transport cross section in eantV internal length2 units
      int indx      = (int)(2.5+fCharge*1.5); // specify that want to get the first tarnsport xsection
      double t1xsec = (gPWAXsecTable->GetPWATotalXsecForZet(zet))->GetInterpXsec(efEnergy,elowindx,indx);
      lambtr1 += (theAtomicNumDensityVector[i]*t1xsec);
    }
    if (lambtr1>0.0) {
      lambtr1 =1./lambtr1;
    }
  } else {
    // Get SCREENING FROM MOLIER
    // below 1 keV it can give bananas so prevet (1 keV)
    if (efEnergy<geant::keV) {
      efEnergy = geant::keV;
    }
    // total mometum square in geantV internal energy2 units which is GeV2
    double pt2     = efEnergy*(efEnergy+2.0*geant::kElectronMassC2);
    double beta2   = pt2/(pt2+geant::kElectronMassC2*geant::kElectronMassC2);
    int    matindx = theMaterial->GetIndex();
    double bc      = gGSTable->GetMoliereBc(matindx);
    double scra    = gGSTable->GetMoliereXc2(matindx)/(4.0*pt2*bc);
    // total elastic mean free path in geantV internal lenght units
    double lambel  = beta2/bc;//*(1.+fScrA);  //NOTE: we should apply correction to the lambtr1 as well !
    double g1      = 2.0*scra*((1.0+scra)*std::log(1.0/scra+1.0)-1.0);
    lambtr1        = lambel/g1;
  }
  return lambtr1;
}


// NOTE: we changed this !
double GSMSCModel::RandomizeTrueStepLength(Geant::GeantTaskData *td, double tlimit) {
  double tempTLimit = tlimit;
  do {
    tempTLimit = td->fRndm->Gauss(tlimit,0.1*tlimit);
  } while ((tempTLimit<0.) || (tempTLimit>2.*tlimit));
  return tempTLimit;
}


}     // namespace geantphysics
