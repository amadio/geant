#ifndef COMPTON_CROSS_SECTION_H
#define COMPTON_CROSS_SECTION_H
//  physics processes for Gammas
// - Compton Scattering

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

// class TGeoVolume;
// class GeantTrack;

class TGeoElement; 
class TGeoMaterial; 
// class ParticleType;

class ComptonCrossSection  // : TNamed()
{
// methods
 public:
   Double_t ComputeCrossSectionPerAtom( double kineticEnergy, double Z); // , double A); 
   Double_t CrossSectionForMaterial(const TGeoMaterial& tMaterial,
                                    Double_t kinEnergyT); 

   const TGeoElement* SelectRandomAtom(const TGeoMaterial& tMaterial,
                                       Double_t kinEnergyT ); 
                                   //  Double_t tcut, Double_t tmax);

   ComptonCrossSection(); 
   ~ComptonCrossSection(); 

// 
 private: 
   enum        { MaxElements= 100 }; 
   Double_t    fNumAtomsPerVolume[MaxElements];  

   Int_t       fNoSec;  
   Double_t    fPartialXsec[MaxElements];   // First trial ( 1 track only - TOFIX ) 
   // std::vector<Double_t> fPartialXsec;   // Size must be >= MaxElements * Max(NumTracks)

   ClassDef( ComptonCrossSection, 1 ) 
};

#endif
