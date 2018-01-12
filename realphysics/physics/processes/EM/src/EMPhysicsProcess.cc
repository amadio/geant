
#include "EMPhysicsProcess.h"

#include "MaterialCuts.h"
#include "Particle.h"
#include "PhysicsParameters.h"
#include "EMModel.h"
#include "EMModelManager.h"
#include "ELossTableRegister.h"
#include "ELossTableManager.h"

#include "LightTrack.h"

#include <iostream>

namespace geantphysics {

EMPhysicsProcess::EMPhysicsProcess(const std::string &name) : PhysicsProcess(name), fModelManager(nullptr),
  fFinalRange(-1.), fDRoverRange(-1.), fLowestKineticEnergy(-1.), fLinearEnergyLossLimit(1.) {
  // set process type
  SetType(ProcessType::kElectromagnetic);
  // create its EMModelManager member
  fModelManager = new EMModelManager();
}


EMPhysicsProcess::~EMPhysicsProcess() {
  // Will delete all EMModel-s as well that has been added to this EMPhysicsProcess (or more exactly to the
  // fModelManager EMModelManager member of this EMPhysicsProcess.
  delete fModelManager;
}


void EMPhysicsProcess::Initialize() {
  // call the base PhysicsProcess Initialize method: it will be checked there if the process is assigned only to
  // particles that the process is allowed to.
  PhysicsProcess::Initialize();
  // init the model manager that will init the models as well togeter with setting reagions where they active
  // the default active regions are determined by the active regions of the process that is given as parameter
  fModelManager->Initialise(GetListActiveRegions(), GetName(), GetListParticlesAssignedTo());
  if (fModelManager==0) {
    std::cerr<<" **** ERROR: EMPhysicsProcess::Initialize() \n"
             <<"      No any models for Process = "<<GetName()
             <<std::endl;
    exit(-1);
  }
  // if it is an energy loss process register it for the list of particles that the process is assigned to into the
  // energy loss table manager (particle, ,this)
  if (ProcessType::kEnergyLoss==GetType()) {
    for (unsigned long i=0; i<GetListParticlesAssignedTo().size(); ++i) {
      ELossTableRegister::Instance().RegisterEnergyLossProcessForParticle(GetListParticlesAssignedTo()[i]->GetInternalCode(),this);
    }
    // set continuous step limit related paramaters from the PhysicsParameters that active in the same regions as this
    std::vector<bool> &active = GetListActiveRegions();
    int regionIndx = 0;
    for (; regionIndx<int(active.size()) && !active[regionIndx]; ++regionIndx) {}
    // TODO: check: if the processes is inactive everywhere!
    const PhysicsParameters *physPar = PhysicsParameters::GetPhysicsParametersForRegion(regionIndx);
    fDRoverRange = physPar->GetDRoverRange();
    fFinalRange  = physPar->GetFinalRange();
    // set the fLowestKineticEnergy
    fLowestKineticEnergy  = physPar->GetLowestElectronTrackingEnergy();
    // set the linear energy loss limit fraction parameter
    fLinearEnergyLossLimit = physPar->GetLinearEnergyLossLimit();
    //
    // print-out
    std::cerr<<"        EMPhysicsProcess Name = " << GetName()
             << "  is an EnergyLoss process! Registered for particles: ";
    for (unsigned long i=0; i<GetListParticlesAssignedTo().size(); ++i) {
      std::cerr<< GetListParticlesAssignedTo()[i]->GetName();
      if (i<GetListParticlesAssignedTo().size()-1)
        std::cout<<", ";
      else
        std::cout<<std::endl;
    }
  }
  std::cerr<<"        EMPhysicsProcess Name = " << GetName() << "  is initialized! "<< std::endl;
}


double EMPhysicsProcess::ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle *particle,
                                     bool istotal) {
  double dedx = 0.0;
  if (GetType()==ProcessType::kEnergyLoss) {
    // loop over the EMModel-s that are active in the region that the current MaterialCuts belongs to;
    // ask them to provide their dedx contribution;
    // use smoothing between models
    const std::vector<EMModel*> models    = fModelManager->GetModelListInRegion(matcut->GetRegionIndex());
    int   numModels                       = models.size();
//      std::cerr<<"  numModels = "<< numModels<<std::endl;
    if (numModels==0)
      return dedx;
    if (numModels==1) {
      dedx = models[0]->ComputeDEDX(matcut, kinenergy, particle, istotal);
    } else {
      // since we are here we have at least 2 models
      int i=numModels-1;
      for (; i>0 && kinenergy<models[i]->GetLowEnergyUsageLimit(); --i) {}
      double delta = 0.0;
      if (i>0) {
        double lowEnergyLimit = models[i]->GetLowEnergyUsageLimit();
        double dedx1 = models[i-1]->ComputeDEDX(matcut, lowEnergyLimit, particle, istotal);
        double dedx2 = models[i]->ComputeDEDX(matcut, lowEnergyLimit, particle, istotal);
       //std::cout<< "====== " <<kinenergy<<" "<<lowEnergyLimit << dedx1/(geant::MeV/geant::mm) <<"  "<<dedx2/(geant::MeV/geant::mm) <<std::endl;
        if (dedx2>0.0) {
          delta = (dedx1/dedx2-1.0)*lowEnergyLimit;
        }
//          std::cout<< "====== " <<kinenergy<<" "<<delta<<" "<<lowEnergyLimit << dedx1/(geant::MeV/geant::mm) <<"  "<<dedx2/(geant::MeV/geant::mm) <<std::endl;
      }
      dedx  = models[i]->ComputeDEDX(matcut, kinenergy, particle, istotal);
//        std::cerr<< "  dedx = " << dedx/(geant::MeV/geant::mm)<< " delta = "<<delta<<std::endl;
      dedx *= (1.0+delta/kinenergy);
    }
    if (dedx<0.0)
      dedx = 0.0;
  }
  return dedx;
}


double EMPhysicsProcess::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                    const Particle *particle, double /*mass*/) const {
  double xsec = 0.0;
  // loop over the EMModel-s that are active in the region that the current MaterialCuts belongs to;
  // ask them to provide their xsec contribution;
  // use smoothing between models
  const std::vector<EMModel*> models    = fModelManager->GetModelListInRegion(matcut->GetRegionIndex());
  int   numModels                       = models.size();
//      std::cerr<<"  numModels = "<< numModels<<std::endl;
  if (numModels==0)
    return xsec;
  if (numModels==1) {
    xsec = models[0]->ComputeMacroscopicXSection(matcut, kinenergy, particle);
  } else {
    // since we are here we have at least 2 models
    int i=numModels-1;
    for (; i>0 && kinenergy<models[i]->GetLowEnergyUsageLimit(); --i) {}
    double delta = 0.0;
    if (i>0) {
      double lowEnergyLimit = models[i]->GetLowEnergyUsageLimit();
      double xsec1 = models[i-1]->ComputeMacroscopicXSection(matcut, lowEnergyLimit, particle);
      double xsec2 = models[i]->ComputeMacroscopicXSection(matcut, lowEnergyLimit, particle);
     //std::cout<< "====== " <<kinenergy<<" "<<lowEnergyLimit << dedx1/(geant::MeV/geant::mm) <<"  "<<dedx2/(geant::MeV/geant::mm) <<std::endl;
      if (xsec2>0.0) {
        delta = (xsec1/xsec2-1.0)*lowEnergyLimit;
      }
//          std::cout<< "====== " <<kinenergy<<" "<<delta<<" "<<lowEnergyLimit << dedx1/(geant::MeV/geant::mm) <<"  "<<dedx2/(geant::MeV/geant::mm) <<std::endl;
    }
    xsec  = models[i]->ComputeMacroscopicXSection(matcut, kinenergy, particle);
//        std::cerr<< "  dedx = " << dedx/(geant::MeV/geant::mm)<< " delta = "<<delta<<std::endl;
    xsec *= (1.0+delta/kinenergy);
  }
  if (xsec<0.0)
    xsec = 0.0;

 return xsec;
}


double EMPhysicsProcess::GetMinimumLambdaTableKineticEnergy(const MaterialCuts *matcut, const Particle*) const {
  double emin = 0.;
  EMModel *lowestEnergyModel =  fModelManager->SelectModel(0., matcut->GetRegionIndex());
  if (lowestEnergyModel) {
    emin = lowestEnergyModel->MinimumPrimaryEnergy(matcut, GetParticle());
  }
  return emin;
}

double EMPhysicsProcess::AlongStepLimitationLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData * /*td*/) const {
  double stepLimit = GetAVeryLargeValue();
  // if the process is kEnergyLoss process use the energy loss related data to limit the step
  if(GetType()==ProcessType::kEnergyLoss) {
    void *mcptr = const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetMaterialCutsPtr();
    const MaterialCuts *matCut = static_cast<const MaterialCuts*>(mcptr);
    double ekin                = gtrack->T();
    const Particle     *part   = Particle::GetParticleByInternalCode(gtrack->GVcode());
    double range = ELossTableManager::Instance().GetRestrictedRange(matCut, part, ekin);
    stepLimit    = range;
    if (range>fFinalRange) {
      stepLimit = range*fDRoverRange+fFinalRange*(1.0-fDRoverRange)*(2.0-fFinalRange/range);
    }
  }
  return stepLimit;
}


//
// Hereve go: the range needs to be re-computed here. It was obtaind at the along step limit but we cannot store!!!
//            so we need to get it again!!!
//
// We need to check on the caller side if kinetic energy become zero after this call and we need to set the track
// status track.SetTrackStatus(??); should be set to stopped but alive i.e. we should check here is the partcile
// has anything to do at rest
int EMPhysicsProcess::AlongStepDoIt(LightTrack &track, Geant::GeantTaskData * /*td*/) {
  int numSecondaries = 0;
//  if (!IsApplicable(track)) {
//    return numSecondaries;
//  }
  if (GetType()==ProcessType::kEnergyLoss) {
    double stepLength = track.GetStepLength(); // this is the current true step length
    // zero step length => nothing to do
    if (stepLength<=0.0) {
      return numSecondaries;
    }
    // check if the partcile is stopping:
    //  - the step length >= range of the particle
    //  - the pre-step point kinetic energy of the partcile is below the common charge partcile trakin cut
    const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
    const Particle     *part   = Particle::GetParticleByInternalCode(track.GetGVcode());
    double              ekin   = track.GetKinE();
    double range = ELossTableManager::Instance().GetRestrictedRange(matCut, part, ekin);
    if (stepLength>=range || ekin<=fLowestKineticEnergy) {
      // the partcile kinetic energy goes to energy deposit
      // update primary track
      track.SetKinE(0.0);
      track.SetEnergyDeposit(track.GetEnergyDeposit()+ekin);
//      track.SetTrackStatus(LTrackStatus::kStopButAlive); // kinetic energy will be checked in the caller
      return numSecondaries;
    }
    //
    // need to compute energy loss
    // 1. try linear energy loss approximation
    double eloss = stepLength*ELossTableManager::Instance().GetRestrictedDEDX(matCut, part, ekin);
    // 2. check if linear approximation is fine i.e. compare to the allowed energy loss fraction and
    //    use integral value if needed
    if (eloss>ekin*fLinearEnergyLossLimit) {
      double postStepRange = range-stepLength;
      eloss = ekin-ELossTableManager::Instance().GetEnergyForRestrictedRange(matCut, part, postStepRange);
    }
    //
    // check final kinetic energy:
    // - if the final energy of the partcile goes below the common tracking cut need to stop the partcile
    double finalEkin = ekin-eloss;
    if (finalEkin<=fLowestKineticEnergy) {
      eloss     = ekin; // everything goes into energy deposit
      finalEkin = 0.0;
//      track.SetTrackStatus(LTrackStatus::kStopButAlive);  kinetic energy will be checked in the caller
    }
    eloss = std::max(eloss, 0.0);
//    std::cout<<" eloss = "<<eloss<<std::endl;
    track.SetKinE(finalEkin);
    track.SetEnergyDeposit(track.GetEnergyDeposit()+eloss);
  }
  return numSecondaries;
}

//
//
//   // set number of interaction left to -1.0 to indicate that new sampling is needed in the next step discrete step limit
//      track.SetNumOfInteractionLegthLeft(-1.0);
// in case of integral approach we need to sample if interaction happens or not; the total macroscopic xsection that
// was used to sample or update the current discrete interaction point is in the track
// if (/*fIsUseIntegralApproach*/) {
//   // get the current total macroscopic cross section i.e. for the post step energy one
//   const Particle *part = Particle::GetParticleByInternalCode(track.GetGVcode());
//   double curInvLambda  =
// }
// IN case of intergral approach it is called only if the discrete interaction indeed happens
int EMPhysicsProcess::PostStepDoIt(LightTrack &track , Geant::GeantTaskData *td) {
  int numSecondaries = 0;
  // for kEnergyLoss processes: check if the particle energy is below the common tracking cut and do nothing if this is
  //                            the case
  double ekin   = track.GetKinE();
  if (GetType()==ProcessType::kEnergyLoss && ekin<=fLowestKineticEnergy) {
    return numSecondaries;
  }
  // ask the model manager to select one model: if the current kinetic energy is below the lowest enedgy model it will
  // return with the lowest energy model => in the SampleSecondaries method of the model we will return
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
  EMModel* model             = fModelManager->SelectModel(ekin, matCut->GetRegionIndex());
  if (model) {
    numSecondaries = model->SampleSecondaries(track, td);
  }
  return numSecondaries;
}


int EMPhysicsProcess::AddModel(EMModel *model) {
  int indx = -1;
  if (!model) {
    std::string partNames = "\n";
    for (unsigned long p=0; p<GetListParticlesAssignedTo().size(); ++p)
      partNames += GetListParticlesAssignedTo()[p]->GetName() + "  ";
    partNames += "\n";
    std::cerr<<" *** WARNING: EMPhysicsProcess::AddModel() \n"
             <<"  Attempt to add EMModel = nullptr to process with Name = " << GetName() << " that is assigned to particles:"
             << partNames << " is ignored."
             << std::endl;
  } else {
    indx = fModelManager->AddModel(model);
  }
  return indx;
}


std::ostream& operator<<(std::ostream& flux,  EMPhysicsProcess* proc) {
  flux << "  Process name = "<<proc->GetName() << " active in regions:  ";
  for (unsigned long i=0; i<proc->GetListActiveRegions().size(); ++i) {
    flux << "\n  R [" << i <<"] = ? ==> " << proc->GetListActiveRegions()[i] << "\n   Models for this region:";
    std::vector<EMModel*> models = proc->fModelManager->GetModelListInRegion(i);
    for (unsigned long j=0; j<models.size(); ++j) {
      flux<<"\n    --> Name = "<<models[j]->GetName() << "  Elow = " <<models[j]->GetLowEnergyUsageLimit()
          << "  Ehigh = " <<models[j]->GetHighEnergyUsageLimit();
      std::vector<bool> actv = models[j]->GetListActiveRegions();
      flux << "\n      ====  model is active in regions: ";
      for (unsigned long k=0; k<actv.size(); ++k)
        flux << "  R[" << k << "] = ? ==> " <<actv[k];

    }
  }
  return flux;
}


} // namespace geantphysics
