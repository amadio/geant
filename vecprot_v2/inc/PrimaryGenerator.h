//===--- PrimaryGenerator.h - Geant-V ---------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file PrimaryGenerator.h
 * @brief Implementation of primary generators for Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANTV_PrimaryGenerator_h
#define GEANTV_PrimaryGenerator_h

class TParticlePDG;
class GeantTrack;

/**
 * @brief Class of primary generators
 */
class PrimaryGenerator : public TNamed {

public:

 /**
  * @brief Pure virtual function of initialization of primary generator
  * @details Set one GeantTrack primary track properties
  */
  virtual void InitPrimaryGenerator() = 0;

  /** @brief  Pure virtual function that produce next event */
  virtual Int_t NextEvent() = 0;

 /**
  * @brief Pure virtual function that returns track
  * 
  * @param n ?????
  * @param gtrack track
  */
  virtual void GetTrack(Int_t n, GeantTrack &gtrack) = 0;
};

#endif
