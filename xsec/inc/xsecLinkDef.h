/* @(#)root/base:$Id: LinkDef2.h 42250 2011-11-25 17:30:31Z pcanal $ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link C++ class TEFstate+;
#pragma link C++ class TEXsec+;
#pragma link C++ class TFinState+;
#pragma link C++ class TMXsec+;
#pragma link C++ class TPDecay+;
#pragma link C++ class TPFstate-;
#pragma link C++ class TPXsec-;
#pragma link C++ class TPartIndex-;
#pragma link C++ class TTabPhysMgr+;
#pragma link C++ class TTabPhysProcess+;
#pragma link C++ class TPrimaryGenerator+;
#pragma link C++ class GunGenerator+;

#ifdef HEPMC
#pragma link C++ class HepMCGenerator+;
#endif

//#pragma link C++ class vecgeom::Vector<TEXsec*>+;

#endif
