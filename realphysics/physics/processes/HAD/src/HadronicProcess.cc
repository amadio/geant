#include "HadronicProcess.h"
#include "LightTrack.h"
#include "HadronicCrossSectionStore.h"
#include "HadronicFinalStateModelStore.h"
#include "Isotope.h"

using namespace geantphysics;

//-----------------------------------
// HadronicProcess non-inline methods
//-----------------------------------
 
HadronicProcess::HadronicProcess() : PhysicsProcess(""), fType( HadronicProcessType::kNotDefined ), 
  fXsecStore( nullptr ), fModelStore( nullptr )
{}


HadronicProcess::HadronicProcess( const std::string &name, const std::vector< int > &particlecodevec, 
                                  const HadronicProcessType type, const bool isatrest,
                                  HadronicCrossSectionStore* xsecstore, HadronicFinalStateModelStore* modelstore ) :
  PhysicsProcess( true, false, isatrest, ForcedCondition::kNotForced, ProcessType::kHadronic, name ),
  fType( type ), fXsecStore( xsecstore ), fModelStore( modelstore )
{
  SetParticleCodeVec( particlecodevec );
}


HadronicProcess::~HadronicProcess() {}


bool HadronicProcess::IsApplicable( const int particlecode ) const {
  bool isOK = false;
  for ( int i = 0; i < fParticleCodeVec.size(); i++ ) {
    if ( fParticleCodeVec[i] == particlecode ) {
      isOK = true;
      break;
    } 
  }
  return isOK;
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
  Material* material = nullptr;
  //***LOOKHERE*** TO-DO : from the indexMaterialCutCouple get the material
  if ( fXsecStore ) {
    std::pair< int, int > pairZandN = fXsecStore->SampleTarget( particleCode, eKin, track.GetMass(), material );
    track.SetTargetZ( pairZandN.first );
    track.SetTargetN( pairZandN.second );
    targetIsotope = Isotope::GetIsotope( pairZandN.first, pairZandN.second );
  }
  return targetIsotope;
}


void HadronicProcess::PostStepDoIt( LightTrack &track, std::vector< LightTrack* > &output ) const {

  // This method does the Lorentz boost of the primary from the Lab to the center-of-mass frame,
  // and the 3D spatial rotation to bring the primary direction from the initial arbitrary one to the z-axis.
  // The first argument is kept constant; the other three are the output of method.
  // QUESTION: IS IT A GOOD IDEA TO ASSUME THIS TRANSFORMATION FOR ALL HADRONIC PROCESSES, INCLUDING ELASTIC?
  //BoostFromLabToCmsAndRotateToMakePrimaryAlongZ( track, transformedTrack, boost, rotation );
  
  //Isotope* targetIsotope = SampleTarget( track );

  // Call now the hadronic model to get the secondaries:
  //int indexModel = 
  //  GetFinalStateModelStore()->GetIndexChosenFinalStateModel( track.GetGVcode(), track.GetKinE(), targetIsotope );
  //  ( GetFinalStateModelStore()->GetHadronicFinalStateModelVec() )[ indexModel ]->SampleFinalState( track, targetIsotope, output );

  //BoostBackFromCmsToLabAndRotateBackToOriginalPrimaryDirection( output, boost, rotation );

  // This method checks all the conservations - charge, energy, momentum, etc. - between the initial state 
  // (the primary "track" and the target nucleus "target") and the final state (the secondaries "output").
  // This check could be done in the center-of-mass frame, or in the lab-frame.
  // QUESTION: BETTER TO CHECK CONSERVATION IN THE CMS FRAME OR IN THE LAB ?
  //CheckConservations( track, targetIsotope, output );

}

void HadronicProcess::AtRestDoIt( LightTrack &track, std::vector< LightTrack* > &output ) {}
