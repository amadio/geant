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

#include <TNamed.h>
#include <TEXsec.h>

class TMXsec : public TNamed {

public:
   TMXsec();
   TMXsec(Int_t z[], Int_t a[], Float_t w[], Int_t nel, Float_t dens, Bool_t weight=kFALSE);
   ~TMXsec() {}

private:
   static TEXsec **fElements; //! pointer to all loaded elements

   Int_t        *fElems;   // List of elements composing this material
   Int_t         fNElems;  // Number of elements
   Float_t      *fTotXS;   // Total x-sec for this material
   Float_t      *fRelXS;   // Relative x-sec for this material

   ClassDef(TMXsec,1)  //Material X-secs

};


#endif
