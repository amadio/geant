/* @(#)root/base:$Id: LinkDef2.h 42250 2011-11-25 17:30:31Z pcanal $ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link C++ class FastSimProcess+;
#pragma link C++ class Smearer+;
#pragma link C++ class GunGenerator+;

#ifdef HEPMC
#pragma link C++ class HepMCGenerator+;
#endif

#endif
