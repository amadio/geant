// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMXsec
#define ROOT_TMXsec


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMXSec                                                               //
//                                                                      //
// X-section for G5 per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <RTypes.h>
#include <TNamed.h>
class TPXsec;
class TGraph;

class TMXsec : public TNamed {

public:
   TMXsec();
   TMXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen, Int_t np);
   ~TMXsec();
   Bool_t AddPart(Int_t kpart, Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax);
   Bool_t AddPartXS(Int_t kpart, const Float_t xsec[], const Short_t dict[]);
   Bool_t AddPartIon(Int_t kpart, const Float_t dedx[]);
   Bool_t Finalise();
   
   Int_t Mat() const {return fMat;}
   Float_t Emin() const {return fEmin;}
   Float_t Emax() const {return fEmax;}
   Int_t NEbins() const {return fNEbins;}
   Double_t ElDelta() const {return fElDelta;}
   Float_t XS(Int_t pdg, Short_t rcode, Float_t en) const;
   Float_t XSindex(Int_t pindex, Short_t rindex, Float_t en) const;
   TGraph *XSGraph(const char* part, const char *reac, 
		   Float_t emin, Float_t emax, Int_t nbin) const;
   void Dump() const;

private:
   Int_t          fMat;     // Material code Z*10000+A*10+metastable level
   Float_t        fEmin;    // Min en in GeV
   Float_t        fEmax;    // Max en in Gev
   Short_t        fNEbins;  // Number of log steps in energy
   Double_t       fElDelta; // Log energy step
   Int_t          fNpart;   // Number of particles
   TPXsec        *fPXsec;   // [fNpart] Cross section table per particle
   Double_t      *fCuts;    // [fNpart] Just a placeholder for the moment

   ClassDef(TMXsec,1)  //Material X-secs

};


#endif
