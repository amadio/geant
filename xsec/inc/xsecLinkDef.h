/* @(#)root/base:$Id: LinkDef2.h 42250 2011-11-25 17:30:31Z pcanal $ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#ifdef __CLING__
#include <string>
#else
#include "dll_stl/str.h" 	 
#endif

#pragma link C++ class TMXsec+;
#pragma link C++ class TEXsec+;
#pragma link C++ class TPXsec-;
#pragma link C++ class TPFstate-;
#pragma link C++ class TFinState+;
#pragma link C++ class TPartIndex-;

#endif
