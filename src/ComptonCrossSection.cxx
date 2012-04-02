#include "ComptonCrossSection.h"
#include "TGeoMaterial.h" 
// #include "TGeoElement.h" 
#include "TRandom.h" 

#include "physical_constants.h"  

#include <cmath>
#include <iostream> 

using namespace std; 
// using std::cout; 

// ClassImp(ComptonCrossSection)
//
// Implementation of Klein Nishina cross section for Compton process
//   copied and adapted from Geant4 class G4KleinNishinaCompton
//    which was created by Michel Maire and Vladimir Ivantchenko
//
//  Adapted by John Apostolakis,  March 8-20th, 2012 


// Will not use G4 classes during tracking
//  May use helper class that 'talks to G4 classes' to create tables at initialisation
//

ComptonCrossSection::ComptonCrossSection()
{
   fNoSec = MaxElements;             //  Works for 1 track at a time only  TOFIX !
   // fPartialXsec.resize( fNoSec ); //  Needs to be  MaxElements * Max(NumTracks) or resizable !
}


ComptonCrossSection::~ComptonCrossSection()
{
}

//
//   Current conventions:
//           GammaEnergy   is in MeV 
//           cross-section returned will be in the unit of barn
//

Double_t 
ComptonCrossSection::ComputeCrossSectionPerAtom(
                                             Double_t GammaEnergy,
                                             Double_t Z )
{
   // Must deal consistently with units -- these are temporary 
   //                        (TODO: check consistency, create & enforce policy )
   const  Double_t  keV= 0.001 * MeV; 
  // Ensure that Gamma is expressed in MeV -- or revise "keV" definition above
  const  Double_t  barn= 1.0; // Cross-section will be in unit of barn

  Double_t CrossSection = 0.0 ;
  if (( Z < 0.9999 ) || ( GammaEnergy < 0.1*keV ) ) return CrossSection;
  //  if ( GammaEnergy > (100.*GeV/Z) ) return CrossSection;

  static const Double_t a = 20.0 , b = 230.0 , c = 440.0;
  
  static const Double_t
    d1= 2.7965e-1*barn, d2=-1.8300e-1*barn, d3= 6.7527   *barn, d4=-1.9798e+1*barn,
    e1= 1.9756e-5*barn, e2=-1.0205e-2*barn, e3=-7.3913e-2*barn, e4= 2.7079e-2*barn,
    f1=-3.9178e-7*barn, f2= 6.8241e-5*barn, f3= 6.0480e-5*barn, f4= 3.0274e-4*barn;
     
  Double_t  Zsq = Z*Z; 
  Double_t p1Z = Z*(d1 + e1*Z + f1*Zsq), p2Z = Z*(d2 + e2*Z + f2*Zsq),
           p3Z = Z*(d3 + e3*Z + f3*Zsq), p4Z = Z*(d4 + e4*Z + f4*Zsq);

  Double_t T0  = 15.0*keV; 
  if (Z < 1.5) T0 = 40.0*keV; 

  Double_t X   = std::max(GammaEnergy, T0) / electron_mass_c2;
  CrossSection = p1Z*log(1.+2.*X)/X
     + (p2Z + p3Z*X + p4Z*X*X)/(1. + (a + (b+ c*X) *X ) *X);
		
  //  modification for low energy (with special case for Hydrogen)
  if (GammaEnergy < T0) {
    Double_t dT0 = 1.*keV;
    X = (T0+dT0) / electron_mass_c2 ;
    Double_t sigma = p1Z*std::log(1.+2*X)/X
                    + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    Double_t   c1 = -T0*(sigma-CrossSection)/(CrossSection*dT0);             
    Double_t   c2 = 0.150; 
    if (Z > 1.5) c2 = 0.375-0.0556*std::log(Z);   // Can lookup log(Z)
    Double_t    y = std::log(GammaEnergy/T0);
    CrossSection *= exp(-y*(c1+c2*y));          
  }
  //  std::cout << "e= " << GammaEnergy << " Z= " << Z << " cross= " << CrossSection << std::endl;
  return CrossSection;
}


// Adapted from 
// G4double G4VEmModel::CrossSectionPerVolume(const G4Material* material, ..


Double_t 
ComptonCrossSection::CrossSectionForMaterial(const TGeoMaterial& tMaterial,
                                 Double_t kinEnergyT)
{
  Double_t   Xsection=0.0;
  Int_t nelm= 1;
  Int_t numTracks= 1; 

  static const Double_t Avogadro = 6.02214179e+23 * 1000;  // */mole;
  Double_t  density=  tMaterial.GetDensity(); 

   // Loop over the elements of  tMaterial

  //  SetupForMaterial(p, material, ekin);
  // const G4ElementVector* theElementVector = material->GetElementVector();
  // const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  const TGeoMixture* ptMixture=  dynamic_cast<const TGeoMixture*>( &tMaterial ); 
  if( ptMixture ) 
  {
     nelm = tMaterial.GetNelements(); 

     Int_t  numEntries= nelm * numTracks; 
     if(numEntries > fNoSec) {
        // fPartialXsec.resize( numEntries );  // Choice for storing partial sums of cross sections
        // fNoSec = numEntries;
        std::cerr << " ERROR in ComptonCrossSection::CrossSectionForMaterial" 
                  << " Material has too many elements: " << numEntries << "." 
                  << " Maximum expected is " << MaxElements  << std::endl;
        std::cerr << " ABORTING. " << std::endl;
        abort(); 
     }
  }

  for (Int_t i=0; i<nelm; i++) {
     //  Constant for this volume -- and shared between tracks
     fNumAtomsPerVolume[i] = Avogadro* density * ptMixture->GetWmixt()[i] / ptMixture->GetAmixt()[i]; 
  }

  for (Int_t i=0; i<nelm; i++) {
     //  Per track value(s)
     TGeoElement* elementI= tMaterial.GetElement(i); 

     Xsection += 
        fNumAtomsPerVolume[i] *
        ComputeCrossSectionPerAtom( kinEnergyT, elementI->Z() ); // ,emin,emax);
     fPartialXsec[i] = Xsection;
     // Store partial sums, to use with Select Random Atom method (for interaction)
  }
  // Xsection= theG4ComptonProcess->CrossSectionPerVolume( kinEnergyT, pMaterialCutCouple ); 

  return     Xsection;
}

const TGeoElement* 
ComptonCrossSection::SelectRandomAtom(const TGeoMaterial& tMaterial,
                                      Double_t kinEnergyT,
                                      TRandom *rngEngine )
{
  const TGeoElement* tCurrentElement= 0; 

  // Double_t *rndArray = &gPropagator->fDblArray[2*tid*ntracks];   //  --> But make sure it is not used already

  const TGeoMixture* ptMixture=  dynamic_cast<const TGeoMixture*>( &tMaterial ); 
  if( ptMixture ) 
  {
     Int_t  nelm = tMaterial.GetNelements(); 
     Int_t  n = nelm - 1;
     tCurrentElement = ptMixture->GetElement(n); 
     Double_t Xsection= fPartialXsec[n];
     // rngEngine->RndmArray(ntracks, rndArray);  
     if (n > 0) {
        Double_t randUnif=  rngEngine->Rndm(); 
        Double_t x = randUnif * Xsection;  // fPartialXsec[n];
           // CrossSectionPerVolume(tMaterial,kinEnergyT); // ,tcut,tmax);
        for(Int_t i=0; i<n; ++i) {
           if (x <= fPartialXsec[i]) {
              tCurrentElement = ptMixture->GetElement(i); 
              std::cout << " Choosing element " << tCurrentElement->GetName() << std::endl;
              std::cout << "    Rand * Xsec= " << x << std::endl;
              std::cout << "    PartialXsec= " << fPartialXsec[i] << std::endl; 
              std::cout << "    Rand= " << randUnif <<  " Xsec= " << fPartialXsec[n] << std::endl; 
              break;
           }
        }
     }
  }else{
     //  Method did not need to be called in this case ... 
     tCurrentElement= tMaterial.GetElement(); 
  }
  return tCurrentElement;
}

