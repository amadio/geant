#include "PhysicsManagerPerParticle.h"

#include "PhysicsParameters.h"
#include "Material.h"
#include "MaterialCuts.h"

#include "PhysicsProcess.h"
#include "Particle.h"

#include "SystemOfUnits.h"
#include "Spline.h"

#include "LightTrack.h"

#include <cmath>
#include <iostream>


namespace geantphysics {

PhysicsManagerPerParticle::PhysicsManagerPerParticle(const Particle *partcile, const PhysicsParameters *physpars)
  : fParticle(partcile), fPhysicsParameters(physpars), fEnergyGrid(nullptr) {
    fIsLambdaTablesPerMaterial  = true;
    fIsHasTotalLambdaTable      = false;
    fIsHasPerProcessLambdaTable = false;
    fIsHasMSCProcess            = false;
    fIsHasDecayProcess          = false;
    fIsHasElossProcess          = false;

    fNumLambdaTableBins         = -1;
    fMinLambdaTableEnergy       = 0.0;
    fMaxLambdaTableEnergy       = 0.0;
    fLogMinLambdaTableEnergy    = 0.0;
    fEnergyILDelta              = 0.0;
  }

// All process object pointers are stored in the PhysicsProcess and all process object is deleted by calling
// PhysicsProcess::ClearAllProcess() static method.
// So here we need to delete only local things
PhysicsManagerPerParticle::~PhysicsManagerPerParticle() {
  // clear the lambda tables
  ClearLambdaTables();
}


void PhysicsManagerPerParticle::Initialize() {
  // Initialize all processes assigned to the particle
  // Note: processes will initialise their models as well
  std::cerr<<"==================  Particle = "<<fParticle->GetName()<< "  has "<<fProcessVec.size()<<" process."<<std::endl;
  for (unsigned long ip=0; ip<fProcessVec.size(); ++ip) {
    fProcessVec[ip]->Initialize();
  }
  // build lambda tables
  BuildLambdaTables();
}


void PhysicsManagerPerParticle::AddProcess( PhysicsProcess *process ) {
  if (process == nullptr) {
    std::cerr<<"   **** ERROR: PhysicsManagerPerParticle::AddProcess()\n"
             <<"    The process pointer is nullptr! \n";
    exit(-1);
  }
  if (process->GetForcedCondition()==ForcedCondition::kInActivated)
    return;
  // the process pointer will be added to this physics manager per particle
  fProcessVec.push_back(process);
  // push also to sub categories
  bool isForced = (process->GetForcedCondition()!=ForcedCondition::kNotForced);
  if (process->GetIsContinuous()) {
    fAlongStepProcessVec.push_back( process );
  }
  if (process->GetIsDiscrete()) {
    fPostStepCandidateProcessVec.push_back(process);
    if (isForced) {
      fPostStepForcedProcessVec.push_back(process);
    }
  }
  if (process->GetIsAtRest())  {
    fAtRestCandidateProcessVec.push_back(process);
    if (isForced) {
      fAtRestForcedProcessVec.push_back(process);
    }
  }
}


double PhysicsManagerPerParticle::GetInvTotalLambda(const MaterialCuts *matcut, double kinenergy) {
  double totMacrXsec = 0.0;
  if (!fIsHasTotalLambdaTable) {
    if (fIsHasDecayProcess) {
      totMacrXsec = fPostStepCandidateProcessVec[0]->ComputeMacroscopicXSection(matcut, kinenergy, fParticle);
    }
    return totMacrXsec;
  }
  LambdaTable *lambTab = nullptr;
  if (fIsLambdaTablesPerMaterial) {
    lambTab = fLambdaTables[matcut->GetMaterial()->GetIndex()];
  } else {
    lambTab = fLambdaTables[matcut->GetIndex()];
  }
  // return with zero if below or above the min/max lambda table energy or there is no lambTab(should not happen)
  if (kinenergy>=fMinLambdaTableEnergy && kinenergy<=fMaxLambdaTableEnergy && lambTab) {
    double logE     = std::log(kinenergy);
    int    lowEIndx = (int) ((logE-fLogMinLambdaTableEnergy)*fEnergyILDelta);
    if (lowEIndx>=fNumLambdaTableBins-1) --lowEIndx;
    // we might put it under verbose build since
    // protection against very small numerical uncertainties
//      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
//        --lowEIndx;
//      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
//        ++lowEIndx;
//      }
    totMacrXsec = lambTab->fSplineInvTotalLambda->GetValueAt(kinenergy,lowEIndx);
    if (totMacrXsec<0.0) totMacrXsec = 0.0;
  }
  return totMacrXsec;
}

double PhysicsManagerPerParticle::ComputeInvTotalLambdaForStepLimit(const MaterialCuts *matcut, double kinenergy) {
  double totMacrXsec = 0.0;
  // if it does not have along step energy loss => get the normal macroscopic cross section and return
  if (!fIsHasElossProcess) {
    totMacrXsec = GetInvTotalLambda(matcut, kinenergy);
    return totMacrXsec;
  }
  // otherwise:
  if (!fIsHasTotalLambdaTable) {
    if (fIsHasDecayProcess) {
      totMacrXsec = fPostStepCandidateProcessVec[0]->ComputeMacroscopicXSection(matcut, kinenergy, fParticle);
    }
    return totMacrXsec;
  }
  LambdaTable *lambTab = nullptr;
  if (fIsLambdaTablesPerMaterial) {
    lambTab = fLambdaTables[matcut->GetMaterial()->GetIndex()];
  } else {
    lambTab = fLambdaTables[matcut->GetIndex()];
  }
  // return with zero if below or above the min/max lambda table energy or there is no lambTab(should not happen)
  if (kinenergy>=fMinLambdaTableEnergy && kinenergy<=fMaxLambdaTableEnergy && lambTab) {
    double ekin = kinenergy;
    // get the energy of the macroscopic cross section maximum
    double maxOfMacXsecE = lambTab->fLambdaMaxEnergy;
    // if the current kinetic energy is already on the left side of this maximum we provide 1/lambda for that just
    // because we assume that 1/lambda is already decreasing on this side with decareasing energy so we provide an
    // overestimate of 1/lambda
    // if the current kinetic energy is on the right side of this maximum point: more work to give an overestimate:
    if (ekin>maxOfMacXsecE) {
      // compute reduced energy: we assume that 1/lambda is higher at lower energy so we provide an overestimate
      double ekinReduced = 0.8*kinenergy;
      // check if it is still on the right side of the maximum point i.e. if our assumption is fine
      // if not: the reduced energy got to the left side of the maximum so we jumped the maximum so return with the max
      if (ekinReduced<maxOfMacXsecE) {
        totMacrXsec = lambTab->fLambdaMax;
        return totMacrXsec;
      } else {
        // otherwise we are still on the right side of the maximum so provide 1/lambda at this reduced energy
        ekin = ekinReduced;
      }
    }
    // if we did not return earlier then we need to provide 1/lambda at ekin that has been set properly above to ensure
    // that the 1/lambda value at ekin energy will be higher than any 1/lambda values along the step i.e. between the
    // current, pre-step and post-step point energy
    double logE     = std::log(ekin);
    int    lowEIndx = (int) ((logE-fLogMinLambdaTableEnergy)*fEnergyILDelta);
    if (lowEIndx>=fNumLambdaTableBins-1) --lowEIndx;
    // we might put it under verbose build since
    // protection against very small numerical uncertainties
//      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
//        --lowEIndx;
//      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
//        ++lowEIndx;
//      }
    totMacrXsec = lambTab->fSplineInvTotalLambda->GetValueAt(ekin,lowEIndx);
    if (totMacrXsec<0.0) totMacrXsec = 0.0;
  }
  return totMacrXsec;
}


PhysicsProcess* PhysicsManagerPerParticle::SelectDiscreteInteraction(const MaterialCuts *matcut, double kinenergy,
                                                                     double presteplambda, Geant::GeantTaskData *td) {
  PhysicsProcess  *proc = nullptr;
  // if there is no any normal (not decay or msc) discrete interactions then there is no total lambda table but
  // on the top of this, if it has only decay it will be selected (discrete process vector has been orderd such that
  // decay is at index =0)
  if (!fIsHasTotalLambdaTable) {
    if (fIsHasDecayProcess) {
      proc = fPostStepCandidateProcessVec[0];
    }
    return proc;
  }
  // get the LambdaTable
  LambdaTable *lambTab = nullptr;
  if (fIsLambdaTablesPerMaterial) {
    lambTab = fLambdaTables[matcut->GetMaterial()->GetIndex()];
  } else {
    lambTab = fLambdaTables[matcut->GetIndex()];
  }
  // if we are here it means that it has at least one normal discrete interaction. If it has more than one normal or
  // one normal plus decay then the per process lambda table exists. Otherwise it has only one normal discerte process
  // that is at index = 0
  // return with zero if below or above the min/max lambda table energy or there is no lambTab(should not happen)
  if (kinenergy>=fMinLambdaTableEnergy && kinenergy<=fMaxLambdaTableEnergy && lambTab) {
    // for particles that has any energy loss processes we gave an overestimated 1/lambda at the step limit so
    // now first we need to see if any interaction happens or just a delta interaction
    if (fIsHasElossProcess) {
      // get the current i.e. for the post-step kinetic energy 1/lambda and compare to the pre-step i.e. overestimated
      // 1/lambda to see if it is just a delta interaction and return immediately with null process if yes
      double curLambda = GetInvTotalLambda(matcut, kinenergy); // current 1/lambda
      // preStepLambda is lambda and not 1/lambda
      if (curLambda<=0.0 || td->fRndm->uniform()>curLambda*presteplambda) { // delta interaction
        proc = nullptr;
        return proc;
      }
    }
    // select interaction
    int numDProc = lambTab->fInvLambdaPerProcess.size();
    if (numDProc==0) {
      proc = fPostStepCandidateProcessVec[0];
    } else {
      // select one discrete process
      double logE     = std::log(kinenergy);
      int    lowEIndx = (int) ((logE-fLogMinLambdaTableEnergy)*fEnergyILDelta);
      if (lowEIndx>=fNumLambdaTableBins-1) --lowEIndx;
      // we might put it under verbose build since
      // protection against very small numerical uncertainties
//      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
//        --lowEIndx;
//      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
//        ++lowEIndx;
//      }

//
//  This version is faster and assumes that normality is not lost during the interpolation and at least one of
//  the probablity of processes is non-zero(that should be in theory since the total lambda was non zero)
      int    procIndx   = 0;
      double cumulative = lambTab->fSplinesInvLambdaPerProcess[procIndx]->GetValueAt(kinenergy,lowEIndx);
      if (cumulative<0.0) cumulative = 0.0;
      double rndm = td->fRndm->uniform();
      for (; rndm>cumulative && procIndx<numDProc-1;) {
        ++procIndx;
        double val = lambTab->fSplinesInvLambdaPerProcess[procIndx]->GetValueAt(kinenergy,lowEIndx);
        if (val<0.0) val = 0.0;
        cumulative += val;
      }
      proc = fPostStepCandidateProcessVec[procIndx];
    }
  }

//
//  This is the full correct way that takes into account the loss of normality due to interpolation.
//
/*
      std::vector<double> cumulative(numDProc,0.0);
      int    procIndx      = 0;
      double prob          = lambTab->fSplinesInvLambdaPerProcess[procIndx]->GetValueAt(kinenergy,lowEIndx);
      cumulative[procIndx] = prob;
      for (procIndx=1; procIndx<numDProc; ++procIndx) {
        prob = lambTab->fSplinesInvLambdaPerProcess[procIndx]->GetValueAt(kinenergy,lowEIndx);
        cumulative[procIndx] = cumulative[procIndx-1]+prob;
      }
      if (cumulative[numDProc-1]>0.0) { // at least one of them is non-zero
        procIndx = 0;
        double dum = rndm*cumulative[numDProc-1];
        for (; dum>cumulative[procIndx] && procIndx<numDProc-1; ++procIndx){}
        proc = fPostStepCandidateProcessVec[procIndx];
      }
//
*/
  return proc;
}


// NOTE: this method is used only for testing and not used during a normal simulation!
// procindex is the index of the (discrete) process in the fPostStepCandidateProcessVec vector
double PhysicsManagerPerParticle::GetMacroscopicXSectionForProcess(const MaterialCuts *matcut, double kinenergy,
                                                                   size_t procindex) {
  double macXsec = 0.0;
  if (procindex>=fPostStepCandidateProcessVec.size()) {
    return macXsec;
  }
  // if there is no any normal (not decay or msc) discrete interactions then there is no total lambda table
  if (!fIsHasTotalLambdaTable) {
    if (fIsHasDecayProcess && procindex==0) { // only decay
      macXsec = fPostStepCandidateProcessVec[0]->ComputeMacroscopicXSection(matcut, kinenergy, fParticle);
    }
    return macXsec;
  }
  // get the LambdaTable
  LambdaTable *lambTab = nullptr;
  if (fIsLambdaTablesPerMaterial) {
    lambTab = fLambdaTables[matcut->GetMaterial()->GetIndex()];
  } else {
    lambTab = fLambdaTables[matcut->GetIndex()];
  }
  // if we are here it means that it has at least one normal discrete interaction. If it has more than one normal or
  // one normal plus decay then the per process lambda table exists. Otherwise it has only one normal discerte process.
  // Return with zero if below or above the min/max lambda table energy or there is no lambTab(should not happen)
  if (kinenergy>=fMinLambdaTableEnergy && kinenergy<=fMaxLambdaTableEnergy && lambTab) {
    // get the current 1/lambda total: was used for normalisation
    double curLambda = GetInvTotalLambda(matcut, kinenergy); // current 1/lambda
    // select interaction
    size_t numDProc = lambTab->fInvLambdaPerProcess.size();
    if (numDProc==0) { // it has only one so Per-Process lambda table was not built
      macXsec = curLambda;
    } else {
      // select one discrete process
      double logE     = std::log(kinenergy);
      int    lowEIndx = (int) ((logE-fLogMinLambdaTableEnergy)*fEnergyILDelta);
      if (lowEIndx>=fNumLambdaTableBins-1) --lowEIndx;
      // we might put it under verbose build since
      // protection against very small numerical uncertainties
//      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
//        --lowEIndx;
//      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
//        ++lowEIndx;
//      }

//
//  This version is faster and assumes that normality is not lost during the interpolation and at least one of
//  the probablity of processes is non-zero(that should be in theory since the total lambda was non zero)
      double normedMacXsec = lambTab->fSplinesInvLambdaPerProcess[procindex]->GetValueAt(kinenergy,lowEIndx);
      macXsec = normedMacXsec*curLambda;
      if (macXsec<0.0) macXsec = 0.0;
    }
  }

//
//  This is the full correct way that takes into account the loss of normality due to interpolation.
//
/*
      std::vector<double> cumulative(numDProc,0.0);
      int    procIndx      = 0;
      double prob          = lambTab->fSplinesInvLambdaPerProcess[procIndx]->GetValueAt(kinenergy,lowEIndx);
      cumulative[procIndx] = prob;
      for (procIndx=1; procIndx<numDProc; ++procIndx) {
        prob = lambTab->fSplinesInvLambdaPerProcess[procIndx]->GetValueAt(kinenergy,lowEIndx);
        cumulative[procIndx] = cumulative[procIndx-1]+prob;
      }
      if (cumulative[numDProc-1]>0.0) { // at least one of them is non-zero
        procIndx = 0;
        double dum = rndm*cumulative[numDProc-1];
        for (; dum>cumulative[procIndx] && procIndx<numDProc-1; ++procIndx){}
        proc = fPostStepCandidateProcessVec[procIndx];
      }
//
*/
  return macXsec;
}


void PhysicsManagerPerParticle::BuildLambdaTables() {
  ClearLambdaTables();
  CheckForLambdaTable();
  if (fIsHasTotalLambdaTable) {
    InitializeEnergyGrid();
    // depending on if per material or per material cuts initilize the fLambdaTables vector with nullptr
    if (fIsLambdaTablesPerMaterial) {
      fLambdaTables.resize(Material::GetTheMaterialTable().size(),nullptr);
    } else {
      fLambdaTables.resize(MaterialCuts::GetTheMaterialCutsTable().size(),nullptr);
    }
    // loop over materials/material cuts that belongs to regions where this PhysicsManagerPerParticle is active and
    // build the tables
    const std::vector<MaterialCuts*> matCutTable = MaterialCuts::GetTheMaterialCutsTable();
    for (unsigned long i=0; i<matCutTable.size(); ++i) {
      MaterialCuts *matCut = matCutTable[i];
      if (fListActiveRegions[matCut->GetRegionIndex()]) {
        if (fIsLambdaTablesPerMaterial) {
          if (!fLambdaTables[matCut->GetMaterial()->GetIndex()]) {
            BuildOneLambdaTable(matCut,matCut->GetMaterial()->GetIndex());
          }
        } else {
          if (!fLambdaTables[matCut->GetIndex()]) {
            BuildOneLambdaTable(matCut,matCut->GetIndex());
          }
        }
      }
    }
  }
}


void PhysicsManagerPerParticle::CheckForLambdaTable() {
  std::vector<int> indxVect;
  int indxDecay = -1;
  int indxMSC   = -1;
  for (unsigned long i=0; i<fPostStepCandidateProcessVec.size(); ++i) {
    PhysicsProcess *proc = fPostStepCandidateProcessVec[i];
    if (proc->GetType()==ProcessType::kMSC) {
      fIsHasMSCProcess = true;
      indxMSC          = i;
    }
    if (proc->GetType()==ProcessType::kDecay) {
      fIsHasDecayProcess = true;
      indxDecay          = i;
    }
    if (proc->GetType()==ProcessType::kEnergyLoss) {
      fIsHasElossProcess         = true;
      fIsLambdaTablesPerMaterial = false; // i.e. per MaterialCuts
    }
    if (proc->GetType()!=ProcessType::kMSC && proc->GetType()!=ProcessType::kDecay) {
      indxVect.push_back(i);
    }
  }
  // 1.
  if (indxVect.size()>0) {
    fIsHasTotalLambdaTable = true;
    // 2.
    if (indxVect.size()==1) {
      // 2. a.
      if (fIsHasDecayProcess) { // we need to make sure that [0,decay], [1,discrete process] ...
        fIsHasPerProcessLambdaTable = true;
      // 2. b.
      } else { // we need to make sure that the only one discrete process is at [0,discrete process] ...
        fIsHasPerProcessLambdaTable = false;
      }
    // 3.
    } else { // has more than one discrete process so we build per process lambda table
     fIsHasPerProcessLambdaTable = true;
    }
  }
  // 4.
  std::vector<PhysicsProcess*> theCopy(fPostStepCandidateProcessVec.size());
  for (unsigned long i=0; i<fPostStepCandidateProcessVec.size(); ++i) {
    theCopy[i] = fPostStepCandidateProcessVec[i];
  }
  if (fIsHasDecayProcess) {
    fPostStepCandidateProcessVec[0] = theCopy[indxDecay];
    for (unsigned long i=1; i<=indxVect.size(); ++i) {
      fPostStepCandidateProcessVec[i] = theCopy[indxVect[i-1]];
    }
    if (fIsHasMSCProcess) {
      fPostStepCandidateProcessVec[fPostStepCandidateProcessVec.size()-1] = theCopy[indxMSC];
    }
  } else {
    for (unsigned long i=0; i<indxVect.size(); ++i) {
      fPostStepCandidateProcessVec[i] = theCopy[indxVect[i]];
    }
    if (fIsHasMSCProcess) {
      fPostStepCandidateProcessVec[fPostStepCandidateProcessVec.size()-1] = theCopy[indxMSC];
    }
  }
}


void PhysicsManagerPerParticle::BuildOneLambdaTable(const MaterialCuts *matcut, int indx) {
  //allocate space for the LambdaTable struct
  LambdaTable *lambTab      = new LambdaTable();
  lambTab->fNumData         = fNumLambdaTableBins;
  lambTab->fMaterialCut     = matcut;
  lambTab->fMaterial        = matcut->GetMaterial();
  lambTab->fLambdaMax       = -1.0;
  lambTab->fLambdaMaxEnergy = -1.0;
  lambTab->fInvTotalLambda  = new double[lambTab->fNumData]();
  lambTab->fSplineInvTotalLambda = nullptr;
  int numDProc = fPostStepCandidateProcessVec.size();
  if (fIsHasMSCProcess) {
    --numDProc; // no for MSC
  }
  if (fIsHasPerProcessLambdaTable) {
    lambTab->fInvLambdaPerProcess.resize(numDProc,nullptr);
    for (int i=0; i<numDProc; ++i) {
      lambTab->fInvLambdaPerProcess[i] = new double[lambTab->fNumData]();
    }
  }
  for (int i=0; i<lambTab->fNumData; ++i) {
    double ekin = fEnergyGrid[i];
    double tot  = 0.0;
    for (int j=0; j<numDProc; ++j) {
      double xsec = fPostStepCandidateProcessVec[j]->ComputeMacroscopicXSection(matcut, ekin, fParticle);
      // just a final check
      if (xsec<0.0) {
        xsec = 0.0;
      }
      if (fIsHasPerProcessLambdaTable) {
        lambTab->fInvLambdaPerProcess[j][i] = xsec;
      }
      tot += xsec;
    }
    lambTab->fInvTotalLambda[i] = tot;
    if (tot>lambTab->fLambdaMax) {
      lambTab->fLambdaMax       = tot;
      lambTab->fLambdaMaxEnergy = ekin;
    }
  }
  //normalization of per process lambdas
  if (fIsHasPerProcessLambdaTable) {
    for (int i=0; i<lambTab->fNumData; ++i) {
      if (lambTab->fInvTotalLambda[i]>0.0) {
        for (int j=0; j<numDProc; ++j) {
          lambTab->fInvLambdaPerProcess[j][i] /= lambTab->fInvTotalLambda[i];
        }
      }
    }
    // set up per process Splines
    lambTab->fSplinesInvLambdaPerProcess.resize(numDProc,nullptr);
    for (int j=0; j<numDProc; ++j) {
      lambTab->fSplinesInvLambdaPerProcess[j] = new Spline(fEnergyGrid, lambTab->fInvLambdaPerProcess[j], lambTab->fNumData);
    }
  }
  lambTab->fSplineInvTotalLambda = new Spline(fEnergyGrid, lambTab->fInvTotalLambda, lambTab->fNumData);
  fLambdaTables[indx] = lambTab;
}


void PhysicsManagerPerParticle::InitializeEnergyGrid() {
  fNumLambdaTableBins      = fPhysicsParameters->GetNumLambdaTableBins()+1;
  fMinLambdaTableEnergy    = fPhysicsParameters->GetMinLambdaTableEnergy();
  fMaxLambdaTableEnergy    = fPhysicsParameters->GetMaxLambdaTableEnergy();
  if (!fEnergyGrid) {
    delete [] fEnergyGrid;
    fEnergyGrid = nullptr;
  }
  fEnergyGrid                        = new double[fNumLambdaTableBins]();
  fLogMinLambdaTableEnergy           = std::log(fMinLambdaTableEnergy);
  double delta                       = std::log(fMaxLambdaTableEnergy/fMinLambdaTableEnergy)/(fNumLambdaTableBins-1.0);
  fEnergyILDelta                     = 1.0/delta;
  fEnergyGrid[0]                     = fMinLambdaTableEnergy;
  fEnergyGrid[fNumLambdaTableBins-1] = fMaxLambdaTableEnergy;
  for (int i=1; i<fNumLambdaTableBins-1; ++i) {
    fEnergyGrid[i] = std::exp(fLogMinLambdaTableEnergy+i*delta);
  }
}


void PhysicsManagerPerParticle::ClearLambdaTables() {
  if (fEnergyGrid) {
    delete [] fEnergyGrid;
    fEnergyGrid = nullptr;
  }
  for (unsigned long i=0; i<fLambdaTables.size(); ++i) {
    if (fLambdaTables[i]) {
      LambdaTable *lambTab = fLambdaTables[i];
      if (lambTab->fInvTotalLambda) {
        delete [] lambTab->fInvTotalLambda;
      }
      if (lambTab->fSplineInvTotalLambda) {
        delete lambTab->fSplineInvTotalLambda;
      }
      for (unsigned long j=0; j<lambTab->fInvLambdaPerProcess.size(); ++j) {
        if (lambTab->fInvLambdaPerProcess[j]) {
          delete [] lambTab->fInvLambdaPerProcess[j];
        }
        if (lambTab->fSplinesInvLambdaPerProcess[j]) {
          delete lambTab->fSplinesInvLambdaPerProcess[j];
        }
      }
      lambTab->fInvLambdaPerProcess.clear();
      lambTab->fSplinesInvLambdaPerProcess.clear();
      delete lambTab;
    }
  }
  fLambdaTables.clear();
  fIsLambdaTablesPerMaterial  = true;
  fIsHasTotalLambdaTable      = false;
  fIsHasPerProcessLambdaTable = false;
  fIsHasMSCProcess            = false;
  fIsHasDecayProcess          = false;
  fIsHasElossProcess          = false;
}


void PhysicsManagerPerParticle::PrepareForRun() {
  // keep only one kEnergyLoss process pointer in the fAlongStepProcessVec if there are more than 1
  if (fIsHasElossProcess) {
    bool isElossProcessFound = false;
    for (unsigned long i=0; i<fAlongStepProcessVec.size(); ++i) {
      if (fAlongStepProcessVec[i]->GetType()==ProcessType::kEnergyLoss) {
        if (!isElossProcessFound) {
          isElossProcessFound = true;
        } else { // remove the process pointer if other kEnergyLoss process had been found before
          fAlongStepProcessVec.erase(fAlongStepProcessVec.begin()+i);
        }
      }
    }
  }
}


void PhysicsManagerPerParticle::ComputeIntLen(LightTrack &track, Geant::GeantTaskData *td) {
  // set total mfp to 1.0 to handle properly particles that has no interaction when
  // fNintLen is updated in the WorkloadManager i.e. step-length/mfp
  track.SetTotalMFP(1.0);// fIntLenV[i] = 1.0;
  //
  // 1. First the continuous step limit if any
  double cStepLimit = PhysicsProcess::GetAVeryLargeValue();
  // if it has MSC than go to last but one
  unsigned long maxContIndx = fAlongStepProcessVec.size();
  if (fIsHasMSCProcess) {
    --maxContIndx;
  }
  for (unsigned long i=0; i<maxContIndx; ++i) {
    //if (fAlongStepProcessVec[i]->IsApplicable(track)) {
      double limit = fAlongStepProcessVec[i]->AlongStepLimitationLength(track);
      if (limit<cStepLimit) {
        cStepLimit = limit;
      }
    //}
  }
  track.SetStepLength(cStepLimit);
  // set fEindeV[] = -1.0 to indicate that continuous step limit has happen; will be updated below if not
  track.SetTargetZ(-1.0);
  //
  // 2. Discrete step limit if any
  //  - check if it has total lambda table i.e. inverse of the total mac. cros. section
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
  double ekin                = track.GetKinE();
  if (fIsHasTotalLambdaTable) { // use that lambda to either update or resample
//
// Without taking into account energy loss along the step
//    double curLambda = GetInvTotalLambda(matCut, ekin); // INTEGRAL APPROACH if it has kENergyLoss process!!!
// With taking into account energy loss along the step if any
    double curLambda = ComputeInvTotalLambdaForStepLimit(matCut, ekin);

    if (curLambda>0.0) {
      curLambda = 1.0/curLambda;
    } else {
      curLambda = PhysicsProcess::GetAVeryLargeValue();
    }
    // check if we need to sample new num.-of-int.-length-left: it is updated in the WorkloadManager after propagation
    if (track.GetNumOfInteractionLegthLeft()<=0.0) {
      double rndm = td->fRndm->uniform(); // use vecgeom RNG to get uniform random number
      track.SetNumOfInteractionLegthLeft(-std::log(rndm));
    }
    track.SetTotalMFP(curLambda);   // save the inverse total mfp; to be used for the update in WorkloadManager
    //update the step length => length = lambda_t * -1. * log(rndm) = lambda_t * number of interaction leght left;
    double limit = curLambda*track.GetNumOfInteractionLegthLeft();
    if (limit<cStepLimit) {
      track.SetStepLength(limit);
      // Flag to indicate that discrete and NOT continous step limit
      track.SetTargetZ(1000);
    }
  } else if (fIsHasDecayProcess) {
    // has no inverse total lambda table so it might has only decay, then it is at index 0
  } // else do nothing
}


int PhysicsManagerPerParticle::AlongStepAction(LightTrack &track, std::vector<LightTrack> &sectracks) {
  int numSecondaries = 0;
  for (unsigned long i=0; i<fAlongStepProcessVec.size(); ++i) {
    numSecondaries += fAlongStepProcessVec[i]->AlongStepDoIt(track, sectracks);
    // check if tarck is still alive
    if (track.GetKinE()<=0.0) { // stopped
      // set track status
      if (fAtRestCandidateProcessVec.size()>0) {
        // invoke AtRestAction;
      } else {
        track.SetTrackStatus(LTrackStatus::kKill);
      }
      return numSecondaries;
    }
  }
  return numSecondaries;
}

int PhysicsManagerPerParticle::PostStepAction(LightTrack &track, std::vector<LightTrack> &sectracks,
                                              Geant::GeantTaskData *td) {
  int numSecondaries = 0;
  // select discrete process
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
  double ekin          = track.GetKinE();
  double preStepLambda = track.GetTotalMFP();
//  double rndm = td->fRndm->uniform(); // use vecgeom RNG to get uniform random number
  PhysicsProcess *proc = SelectDiscreteInteraction(matCut, ekin, preStepLambda, td);
  // should be done like this
//  if (!proc || !proc->IsApplicable(track)) {
//    return numSecondaries;
//  }
  if (!proc) {
    return numSecondaries;
  }
  // invoke the post step action of the selected discrete process
  numSecondaries = proc->PostStepDoIt(track , sectracks, td);
  return numSecondaries;
}


}  // namespace geantphysics
