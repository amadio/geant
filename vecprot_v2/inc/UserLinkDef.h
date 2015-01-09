//===--- UserLinkDef.h - Geant-V --------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file UserLinkDef.h
 * @brief Linkage for classes in user applications
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class MyHit+;
#pragma link C++ class MyApplication+;
