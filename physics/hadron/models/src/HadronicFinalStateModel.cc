#include "HadronicFinalStateModel.h"
#include "Isotope.h"

using namespace geantphysics;

//-------------------------------------------
// HadronicFinalStateModel non-inline methods
//-------------------------------------------

HadronicFinalStateModel::HadronicFinalStateModel() : 
  fName( "" ), fType( HadronicModelType::kNotDefined ), fLowEnergyUsageLimit( 0.0 ), fHighEnergyUsageLimit( 100000.0 ),
  fMinTargetZ( 1 ), fMaxTargetZ( 120 ), fMinTargetN( 1 ), fMaxTargetN( 300 )
{}


HadronicFinalStateModel::
HadronicFinalStateModel( const std::string name, const std::vector< int > &projectilecodevec,
                         const HadronicModelType type, const double minenergy, const double maxenergy,
                         const double mintargetz, const double maxtargetz,
                         const double mintargetn, const double maxtargetn ) :
  fName( name ), fType( type ), fLowEnergyUsageLimit( minenergy ), fHighEnergyUsageLimit( maxenergy ),
  fMinTargetZ( mintargetz ), fMaxTargetZ( maxtargetz ), fMinTargetN( mintargetn ), fMaxTargetN( maxtargetn )
{
  SetProjectileCodeVec( projectilecodevec );
}


HadronicFinalStateModel::HadronicFinalStateModel( const HadronicFinalStateModel &other ) :
  fName( other.fName ), fType( other.fType), fLowEnergyUsageLimit( other.fLowEnergyUsageLimit ), fHighEnergyUsageLimit( other.fHighEnergyUsageLimit ),
  fMinTargetZ( other.fMinTargetZ ), fMaxTargetZ( other.fMaxTargetZ ),
  fMinTargetN( other.fMinTargetN ), fMaxTargetN( other.fMaxTargetN )
{
  SetProjectileCodeVec( other.fProjectileCodeVec );
}


HadronicFinalStateModel& HadronicFinalStateModel::operator=( const HadronicFinalStateModel &other ) {
  if ( this != &other ) {
    SetProjectileCodeVec( other.fProjectileCodeVec );
    fName = other.fName;
    fType = other.fType;
    fLowEnergyUsageLimit = other.fLowEnergyUsageLimit;
    fHighEnergyUsageLimit = other.fHighEnergyUsageLimit;
    fMinTargetZ = other.fMinTargetZ;
    fMaxTargetZ = other.fMaxTargetZ;
    fMinTargetN = other.fMinTargetN;
    fMaxTargetN = other.fMaxTargetN;
  }
  return *this;
}


HadronicFinalStateModel::~HadronicFinalStateModel() {}


void HadronicFinalStateModel::Initialize( /* Not yet defined */ ) {}


bool HadronicFinalStateModel::
IsApplicable( const int projectilecode, const double projectilekineticenergy, const Isotope* targetisotope ) {
  bool isOK = false;
  
  for ( size_t i = 0; i < fProjectileCodeVec.size(); i++ ) {
    if ( fProjectileCodeVec[i] == projectilecode ) {
      isOK = true;
      break;
    } 
  }
  if ( isOK ) {
    if ( projectilekineticenergy < fLowEnergyUsageLimit  ||  projectilekineticenergy > fHighEnergyUsageLimit ) {
      isOK = false;
    } else {
      int targetZ = 0;
      int targetN = 0;
      if ( targetisotope ) {
        targetZ = static_cast< int >( targetisotope->GetZ() );
        targetN = targetisotope->GetN();
      }
      if ( targetZ < fMinTargetZ  ||  targetZ > fMaxTargetZ  ||  targetN < fMinTargetN  ||  targetN > fMaxTargetN ) {
        isOK = false;
      }
    }
  }
  return isOK;
}

