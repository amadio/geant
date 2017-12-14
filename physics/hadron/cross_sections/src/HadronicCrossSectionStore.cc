#include "HadronicCrossSectionStore.h"
#include "HadronicCrossSection.h"
#include "Isotope.h"
#include "Element.h"
#include "Material.h"
#include "MaterialProperties.h"

using namespace geantphysics;

//---------------------------------------------
// HadronicCrossSectionStore non-inline methods
//---------------------------------------------

HadronicCrossSectionStore::HadronicCrossSectionStore() : fName( "" )
{}


HadronicCrossSectionStore::HadronicCrossSectionStore( const std::string name ) : fName( name )
{}


HadronicCrossSectionStore::HadronicCrossSectionStore( const HadronicCrossSectionStore &other ) :
  fName( other.fName )
{
  for ( size_t i = 0; i < other.fHadXsecVec.size(); i++ ) {
    fHadXsecVec.push_back( other.fHadXsecVec[i] );
  }
}


HadronicCrossSectionStore& HadronicCrossSectionStore::operator=( const HadronicCrossSectionStore &other ) {
  if ( this != &other ) {
    fHadXsecVec.clear();
    for ( size_t i = 0; i < other.fHadXsecVec.size(); i++ ) {
      fHadXsecVec.push_back( other.fHadXsecVec[i] );
    }
    fName = other.fName;  
  }
  return *this;
}


HadronicCrossSectionStore::~HadronicCrossSectionStore() {
  // We are assuming here that this class is the owner of the hadronic cross sections, therefore it is in charge
  // of deleting them at the end.
  for ( size_t i = 0; i < fHadXsecVec.size(); i++ ) {
    delete fHadXsecVec[i];
  }
  fHadXsecVec.clear();
}


void HadronicCrossSectionStore::Initialize( /* Not yet defined */ ) {}


void HadronicCrossSectionStore::
RegisterHadronicCrossSection( HadronicCrossSection* ptrhadxsec ) {
  if ( ptrhadxsec ) {
    fHadXsecVec.push_back( ptrhadxsec );
  }
}


int HadronicCrossSectionStore::
GetIndexFirstApplicableXsec( const int projectilecode, const double projectilekineticenergy,
                             const Element* targetelement, const Material* /*targetmaterial*/ ) {
  int index = -1;
  for ( int i = fHadXsecVec.size() - 1; i >= 0 ; i-- ) {
    if ( fHadXsecVec[i]  && 
         fHadXsecVec[i]->IsApplicable( projectilecode, projectilekineticenergy, targetelement->GetZ(), targetelement->GetA()/(geant::g/geant::mole) ) ) {
      index = i;
      break;
    }
  }
  return index;
}


double HadronicCrossSectionStore::
GetIsotopeCrossSection( const int projectilecode, const double projectilekineticenergy, const double projectilemass,
			const Isotope* targetisotope, const Element* targetelement, const Material* targetmaterial) {
  double xsec = -1.0;
  int index = GetIndexFirstApplicableXsec( projectilecode, projectilekineticenergy, targetelement, targetmaterial );
  if ( index >= 0 ) {      
    xsec = fHadXsecVec[index]->GetIsotopeCrossSection( projectilecode, projectilekineticenergy, 
                                                       projectilemass, targetisotope->GetZ(), targetisotope->GetN() );
  }
  return xsec;
}


double HadronicCrossSectionStore::
GetElementCrossSection( const int projectilecode, const double projectilekineticenergy, const double projectilemass,
                        const Element* targetelement, const Material* targetmaterial ) {
  double xsec = -1.0;
  int index = GetIndexFirstApplicableXsec( projectilecode, projectilekineticenergy, targetelement, targetmaterial );

  if ( index >= 0 ) {
    const Vector_t< geantphysics::Isotope* > isotopeVector = targetelement->GetIsotopeVector();
      const double* abundanceIsotopeVector = targetelement->GetRelativeAbundanceVector();
      xsec = 0.0;
      for ( size_t i = 0; i < isotopeVector.size(); i++ ) {

	double isotopeXsec = fHadXsecVec[index]->GetIsotopeCrossSection( projectilecode, projectilekineticenergy, projectilemass,
									 isotopeVector[i]->GetZ(), isotopeVector[i]->GetN());
	if ( isotopeXsec < 0.0 ) {
	  xsec = -1.0;
	  break;
	}
	xsec += abundanceIsotopeVector[i] * isotopeXsec;
    }
  }
  return xsec;
}


double HadronicCrossSectionStore::
GetMacroscopicCrossSection( const int projectilecode, const double projectilekineticenergy, const double projectilemass,
			    const Material* targetmaterial ) {
  double xsec = -1.0;
  if ( targetmaterial ) {
    Vector_t< Element* > elementVector = targetmaterial->GetElementVector();
    const double* numOfAtomsPerVolumeVector = targetmaterial->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
    xsec = 0.0;
    for ( size_t i = 0; i < elementVector.size(); i++ ) {
      double elementXsec = GetElementCrossSection( projectilecode, projectilekineticenergy, projectilemass,
                                                   elementVector[i], targetmaterial );
      if ( elementXsec < 0.0 ) {
        xsec = -1.0;
        continue;
      }
      xsec += elementXsec * numOfAtomsPerVolumeVector[i];
    }
  }
  return xsec;
}


std::pair< int, int > HadronicCrossSectionStore::
SampleTarget( const int projectilecode, const double projectilekineticenergy, const double projectilemass,
	      const Material* targetmaterial ) {
  // On-the-fly sampling of the target isotope.
  // It would be better to replace this with a much faster method that relies on some pre-computed tables built
  // at initialization, during the building of the lambda tables.
  int targetZ = 0;
  int targetN = 0;
  if ( targetmaterial ) {

    // First select the element, i.e. the atomic number Z
    Vector_t< Element* > elementVector = targetmaterial->GetElementVector();
    const double* numOfAtomsPerVolumeVector = targetmaterial->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
    std::vector< double > sumElementXsecVector;
    double xsec = 0.0;
    for ( size_t i = 0; i < elementVector.size(); i++ ) {
      double elementXsec = GetElementCrossSection( projectilecode, projectilekineticenergy, projectilemass,
                                                   elementVector[i], targetmaterial );
      if ( elementXsec < 0.0 ) break;
      xsec += elementXsec * numOfAtomsPerVolumeVector[i];  
      sumElementXsecVector.push_back( xsec );
    }
    if ( xsec > 0.0  &&  sumElementXsecVector.size() == elementVector.size() ) {
      double randomNumber1 = 0.5;  //***LOOKHERE*** TO-BE-REPLACED with a call to a random number generator.
      size_t iEle = 0;
      while ( iEle < elementVector.size()  &&  randomNumber1 > sumElementXsecVector[iEle]/xsec ) {
        iEle++;
      }
      targetZ = static_cast< int >( elementVector[iEle]->GetZ() );     

      // Now select the isotope, i.e. the number of nucleons N
      Vector_t< Isotope* > isotopeVector = elementVector[iEle]->GetIsotopeVector();
      const double* abundanceIsotopeVector = elementVector[iEle]->GetRelativeAbundanceVector();
      std::vector< double > sumIsotopeXsecVector;
      xsec = 0.0;
      for ( size_t i = 0; i < isotopeVector.size(); i++ ) {
        double isotopeXsec = GetIsotopeCrossSection( projectilecode, projectilekineticenergy, projectilemass,
                                                     isotopeVector[i], elementVector[iEle], targetmaterial);
        if ( isotopeXsec < 0.0 ) break;
        xsec += abundanceIsotopeVector[i] * isotopeXsec;
        sumIsotopeXsecVector.push_back( xsec );
      }
      if ( xsec > 0.0  &&  sumIsotopeXsecVector.size() == isotopeVector.size() ) {
        double randomNumber2 = 0.5;  //***LOOKHERE*** TO-BE-REPLACED with a call to a random number generator.
        size_t iIso = 0;
        while ( iIso < isotopeVector.size()  &&  randomNumber2 > sumIsotopeXsecVector[iEle]/xsec ) {
          iIso++;
        }
        targetN = isotopeVector[iIso]->GetN();     
      }
    }
  }
  return std::pair< int, int >( targetZ, targetN );
}

