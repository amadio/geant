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

#pragma link C++ class geant::cxx::TaskData+;
#pragma link C++ class Basket+;
#pragma link C++ class Propagator+;
#pragma link C++ class StdApplication+;
#pragma link C++ class MCTruthMgr+;
#pragma link C++ class geant::GeantConfig+;

#endif
