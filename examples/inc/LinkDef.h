//===--- LinkDef.h - Geant-V ------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file LinkDef.h
 * @brief Linkage for classes  
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class ExN03Application+;
#pragma link C++ class CMSApplication+;
#pragma link C++ class LHCbApplication+;
#pragma link C++ class FastSimApplication+;

#endif
