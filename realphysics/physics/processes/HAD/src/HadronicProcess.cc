#include "HadronicProcess.h"
#include "LightTrack.h"
#include "HadronicCrossSectionStore.h"
#include "HadronicFinalStateModelStore.h"
#include "HadronicFinalStateModel.h"
#include "Isotope.h"
#include "Material.h"
#include "MaterialCuts.h"
#include "MaterialProperties.h"

using namespace geantphysics;

//-----------------------------------
// HadronicProcess non-inline methods
//-----------------------------------
 
HadronicProcess::HadronicProcess() : PhysicsProcess(""), fType( HadronicProcessType::kNotDefined ), 
  fXsecStore( nullptr ), fModelStore( nullptr )
{}


HadronicProcess::HadronicProcess( const std::string &name ) :
  PhysicsProcess( name )
{
  fModelStore = new HadronicFinalStateModelStore();
  fXsecStore = new HadronicCrossSectionStore();
}



HadronicProcess::HadronicProcess( const std::string &name, const std::vector< int > &particlecodevec, 
                                  const HadronicProcessType type, const bool isatrest,
                                  HadronicCrossSectionStore* xsecstore, HadronicFinalStateModelStore* modelstore ) :
  PhysicsProcess( true, false, isatrest, ForcedCondition::kNotForced, ProcessType::kHadronic, name ),
  fType( type ), fXsecStore( xsecstore ), fModelStore( modelstore )
{
  SetParticleCodeVec( particlecodevec );
}


HadronicProcess::~HadronicProcess() {}


bool HadronicProcess::IsApplicable( const LightTrack &/*track*/ ) const {

  /*
  int particlecode = track.GetGVcode();
  
  bool isOK = false;
  for ( size_t i = 0; i < fParticleCodeVec.size(); i++ ) {
    if ( fParticleCodeVec[i] == particlecode ) {
      isOK = true;
      break;
    } 
  }
  */
  return true;
}

double HadronicProcess::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
						   const Particle *particle, double mass) const {
  double xsec = 0.0;
  // compute the macroscopic cross section as the sum of the atomic cross sections weighted by the number of atoms in
  // in unit volume.
  const Material *mat =  matcut->GetMaterial();
  // we will need the element composition of this material
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();
  
  for (int iel=0; iel<numElems; ++iel) {
    xsec += theAtomicNumDensityVector[iel]*GetAtomicCrossSection(particle->GetInternalCode(), kinenergy, mass, theElements[iel], mat);
  }

  return xsec;
}


double HadronicProcess::
GetAtomicCrossSection( const int particlecode, const double particlekineticenergy, const double particlemass,
                       const Element* targetelement, const Material* targetmaterial ) const {
  double xsec = -1.0;
  if ( fXsecStore ) {
    xsec = fXsecStore->GetElementCrossSection( particlecode, particlekineticenergy, particlemass, targetelement, targetmaterial );
  }
  return xsec;
}


Isotope* HadronicProcess::SampleTarget( LightTrack &track ) const {
  Isotope* targetIsotope = nullptr;
  int particleCode = track.GetGVcode();
  double eKin = track.GetKinE();
  int indexMaterialCutCouple = track.GetMaterialCutCoupleIndex();

  const Material* material = MaterialCuts::GetTheMaterialCutsTable()[indexMaterialCutCouple]->GetMaterial();
  
  if ( fXsecStore ) {
    std::pair< int, int > pairZandN = fXsecStore->SampleTarget( particleCode, eKin, track.GetMass(), material );
    track.SetTargetZ( pairZandN.first );
    track.SetTargetN( pairZandN.second );
    targetIsotope = Isotope::GetIsotope( pairZandN.first, pairZandN.second );
  }
  return targetIsotope;
}


int HadronicProcess::PostStepDoIt( LightTrack &track, Geant::GeantTaskData *td) {

  // Comment below not up to date anymore
  //
  // This method does the Lorentz boost of the primary from the Lab to the center-of-mass frame,
  // and the 3D spatial rotation to bring the primary direction from the initial arbitrary one to the z-axis.
  // The first argument is kept constant; the other three are the output of method.
  // QUESTION: IS IT A GOOD IDEA TO ASSUME THIS TRANSFORMATION FOR ALL HADRONIC PROCESSES, INCLUDING ELASTIC?
  //BoostFromLabToCmsAndRotateToMakePrimaryAlongZ( track, transformedTrack, boost, rotation );
  
  Isotope* targetIsotope = SampleTarget( track );

  // Call now the hadronic model to get the secondaries:
  int indexModel = 
    GetFinalStateModelStore()->GetIndexChosenFinalStateModel( track.GetGVcode(), track.GetKinE(), targetIsotope );
  ( GetFinalStateModelStore()->GetHadronicFinalStateModelVec() )[ indexModel ]->SampleFinalState( track, targetIsotope, td);

  //BoostBackFromCmsToLabAndRotateBackToOriginalPrimaryDirection( output, boost, rotation );

  // This method checks all the conservations - charge, energy, momentum, etc. - between the initial state 
  // (the primary "track" and the target nucleus "target") and the final state (the secondaries "output").
  // This check could be done in the center-of-mass frame, or in the lab-frame.
  // QUESTION: BETTER TO CHECK CONSERVATION IN THE CMS FRAME OR IN THE LAB ?
  //CheckConservations( track, targetIsotope, output );
  return 0;
}

void HadronicProcess::AtRestDoIt( LightTrack& /*track*/,  Geant::GeantTaskData * /*td*/ ) {}


void HadronicProcess::AddModel(HadronicFinalStateModel *model) {
  if (!model) {
    std::string partNames = "\n";
    for (unsigned long p=0; p<GetListParticlesAssignedTo().size(); ++p)
      partNames += GetListParticlesAssignedTo()[p]->GetName() + "  ";
    partNames += "\n";
    std::cerr<<" *** WARNING: HadronicProcess::AddModel() \n"
             <<"  Attempt to add HadronicFinalStateModel = nullptr to process with Name = " << GetName() << " that is assigned to particles:"
             << partNames << " is ignored."
             << std::endl;
  } else {
    fModelStore->RegisterHadronicFinalStateModel(model);
  }
  return;
}

void HadronicProcess::AddCrossSection(HadronicCrossSection *xsection) {
  if (!xsection) {
    std::string partNames = "\n";
    for (unsigned long p=0; p<GetListParticlesAssignedTo().size(); ++p)
      partNames += GetListParticlesAssignedTo()[p]->GetName() + "  ";
    partNames += "\n";
    std::cerr<<" *** WARNING: HadronicProcess::AddCrossSection() \n"
             <<"  Attempt to add HadronicCrossSection = nullptr to process with Name = " << GetName() << " that is assigned to particles:"
             << partNames << " is ignored."
             << std::endl;
  } else {
    fXsecStore->RegisterHadronicCrossSection(xsection);
  }
  return;
}
  

