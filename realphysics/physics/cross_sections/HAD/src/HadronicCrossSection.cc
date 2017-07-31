#include "HadronicCrossSection.h"
#include <iostream>
using namespace geantphysics;

//----------------------------------------
// HadronicCrossSection non-inline methods
//----------------------------------------

HadronicCrossSection::HadronicCrossSection() : 
  fName( "" ), fMinEnergy( 0.0 ), fMaxEnergy( 100000.0 ),
  fMinTargetZ( 1 ), fMaxTargetZ( 120 ), 
  fMinTargetN( 1 ), fMaxTargetN( 300 )
{}


HadronicCrossSection::HadronicCrossSection( const std::string name, const std::vector< int > &projectilecodevec,
                                            const double minenergy, const double maxenergy,
                                            const double mintargetz, const double maxtargetz,
                                            const double mintargetn, const double maxtargetn ) :
  fName( name ), fMinEnergy( minenergy ), fMaxEnergy( maxenergy ),
  fMinTargetZ( mintargetz ), fMaxTargetZ( maxtargetz ),
  fMinTargetN( mintargetn ), fMaxTargetN( maxtargetn )
{
  SetProjectileCodeVec( projectilecodevec );
}


HadronicCrossSection::HadronicCrossSection( const HadronicCrossSection &other ) :
  fName( other.fName ), fMinEnergy( other.fMinEnergy ), fMaxEnergy( other.fMaxEnergy ),
  fMinTargetZ( other.fMinTargetZ ), fMaxTargetZ( other.fMaxTargetZ ),
  fMinTargetN( other.fMinTargetN ), fMaxTargetN( other.fMaxTargetN )
{
  SetProjectileCodeVec( other.fProjectileCodeVec );
}


HadronicCrossSection& HadronicCrossSection::operator=( const HadronicCrossSection &other ) {
  if ( this != &other ) {
    SetProjectileCodeVec( other.fProjectileCodeVec );
    fName = other.fName;  
    fMinEnergy = other.fMinEnergy;
    fMaxEnergy = other.fMaxEnergy;
    fMinTargetZ = other.fMinTargetZ;
    fMaxTargetZ = other.fMaxTargetZ;
    fMinTargetN = other.fMinTargetN;
    fMaxTargetN = other.fMaxTargetN;
  }
  return *this;
}


HadronicCrossSection::~HadronicCrossSection() {}


void HadronicCrossSection::Initialize( /* Not yet defined */ ) {}


bool HadronicCrossSection::IsApplicable( const int projectilecode, const double projectilekineticenergy,
                                         const int targetZ, const int targetN ) {
  bool isOK = false;
  for ( unsigned int i = 0; i < fProjectileCodeVec.size(); i++ ) {
    if ( fProjectileCodeVec[i] == projectilecode ) {
      isOK = true;
      break;
    } 
  }

  if ( isOK ) {
    if ( projectilekineticenergy < fMinEnergy  ||  projectilekineticenergy > fMaxEnergy ) {
      isOK = false;
    } else {
      if ( targetZ < fMinTargetZ  ||  targetZ > fMaxTargetZ ) {
        isOK = false;
      } else {
          if ( targetN < fMinTargetN  ||  targetN > fMaxTargetN ) {
            isOK = false;
	  }
      }       
    }
  }
  return isOK;
}

