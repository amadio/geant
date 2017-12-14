#include "HadronicFinalStateModelStore.h"
#include "HadronicFinalStateModel.h"
#include "Isotope.h"

using namespace geantphysics;

//------------------------------------------------
// HadronicFinalStateModelStore non-inline methods
//------------------------------------------------

HadronicFinalStateModelStore::HadronicFinalStateModelStore() : fName( "" )
{}


HadronicFinalStateModelStore::HadronicFinalStateModelStore( const std::string name ) : fName( name )
{}


HadronicFinalStateModelStore::HadronicFinalStateModelStore( const HadronicFinalStateModelStore &other ) :
  fName( other.fName )
{
  for ( size_t i = 0; i < other.fHadFsVec.size(); i++ ) {
    fHadFsVec.push_back( other.fHadFsVec[i] );
  }
}


HadronicFinalStateModelStore& HadronicFinalStateModelStore::operator=( const HadronicFinalStateModelStore &other ) {
  if ( this != &other ) {
    fHadFsVec.clear();
    for ( size_t i = 0; i < other.fHadFsVec.size(); i++ ) {
      fHadFsVec.push_back( other.fHadFsVec[i] );
    }
    fName = other.fName;  
  }
  return *this;
}


HadronicFinalStateModelStore::~HadronicFinalStateModelStore() {
  // We are assuming here that this class is the owner of the hadronic final-state models, therefore it is in charge
  // of deleting them at the end.
  for ( size_t i = 0; i < fHadFsVec.size(); i++ ) {
    delete fHadFsVec[i];
  }
  fHadFsVec.clear();
}


void HadronicFinalStateModelStore::Initialize( /* Not yet defined */ ) {}


void HadronicFinalStateModelStore::
RegisterHadronicFinalStateModel( HadronicFinalStateModel* ptrhadfs ) {
  if ( ptrhadfs ) {
    fHadFsVec.push_back( ptrhadfs );
  }
}


int HadronicFinalStateModelStore::
GetIndexChosenFinalStateModel( const int projectilecode, const double projectilekineticenergy,
                               const Isotope* targetisotope ) const {
  int index = -1;
  std::vector< int > indexApplicableModelVec;

  for ( size_t i = 0; i < fHadFsVec.size(); i++ ) {
    if ( fHadFsVec[i]  &&  fHadFsVec[i]->IsApplicable( projectilecode, projectilekineticenergy,
						       targetisotope )) {
      indexApplicableModelVec.push_back( i );
    }
  }
  
  if ( indexApplicableModelVec.size() == 1 ) {
    index = indexApplicableModelVec[0];
  } else if ( indexApplicableModelVec.size() == 2 ) {
    // The "first" index corresponds to the model with the lowest minimal energy
    int first = indexApplicableModelVec[0];
    int second = indexApplicableModelVec[1];
    if ( fHadFsVec[ first ]->GetLowEnergyUsageLimit() > fHadFsVec[ second ]->GetLowEnergyUsageLimit() ) {
      first = second;
      second = indexApplicableModelVec[0];
    }
    if ( fHadFsVec[ first ]->GetHighEnergyUsageLimit() >= fHadFsVec[ second ]->GetHighEnergyUsageLimit() ) {
      std::cerr << "HadronicFinalStateModelStore::GetIndexChosenFinalState : projectilecode=" << projectilecode
                << " ; projectilekineticenergy=" << projectilekineticenergy << " GeV; targetisotope: Z=" 
                << targetisotope->GetZ() << " N=" << targetisotope->GetN() << std::endl
                << "\t NOT allowed full overlapping between two models: " << fHadFsVec[ first ]->GetName()
                << " , "  << fHadFsVec[ second ]->GetName() << std::endl;
    } else {  // All if fine: first model applicable to lower energies than the second model
      // Select one of the two models with probability depending linearly on the projectilekineticenergy
      double probFirst = ( fHadFsVec[ first ]->GetHighEnergyUsageLimit() - projectilekineticenergy ) / 
                         ( fHadFsVec[ first ]->GetHighEnergyUsageLimit() - fHadFsVec[ second ]->GetLowEnergyUsageLimit() );
      double randomNumber = 0.5;  //***LOOKHERE*** TO-BE-REPLACED with a call to a random number generator.
      if ( randomNumber < probFirst ) {
        index = first;
      } else {
        index = second;
      }
    }
  } else {
    std::cerr << "HadronicFinalStateModelStore::GetIndexChosenFinalState : projectilecode=" << projectilecode
              << " ; projectilekineticenergy=" << projectilekineticenergy << " GeV; targetisotope: Z=" 
              << targetisotope->GetZ() << " N=" << targetisotope->GetN() << std::endl
              << "\t wrong number of applicable final-state models: " << indexApplicableModelVec.size() << std::endl;
  }
  return index;
}

