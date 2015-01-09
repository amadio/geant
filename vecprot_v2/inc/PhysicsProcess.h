//===--- PhysicsProcess.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file PhysicsProcess.h
 * @brief Definition of physical processes used in Geant-V prototype
 * @details Toy physics processes for our propagator prototype.
 * Currently including:
 * 1. Single scattering as a discrete process;
 * 2. Energy loss as continuous process;
 * 3. Generic interaction as discrete process, producing secondaries.
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PHYSICSPROCESS
#define GEANT_PHYSICSPROCESS
#include "Geant/Config.h"

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TGeoMaterial;
class GeantTrack;
class GeantTrack_v;
class TGenPhaseSpace;

#include "TMutex.h"

/**
 * @brief Class describing physics processes
 */
class PhysicsProcess : public TNamed {
public:

  /**
   * @enum EProcessType
   * @brief Enum ProcessType definition
   */
  enum EProcessType { kDiscrete = BIT(14), kContinuous = BIT(15) };

public:

  /**
   * @brief PhysicsProcess constructor 
   */
  PhysicsProcess() : TNamed() {}

  /**
   * @brief PhysicsProcess parametrized constructor
   * 
   * @param name Name of physics process
   */
  PhysicsProcess(const char *name) : TNamed(name, "") {}

  /** @brief PhysicsProcess destructor */
  virtual ~PhysicsProcess() {}

  /**
   * @brief Function that check type
   * 
   * @param type EProcessType type
   * @return Boolean value -> (Bool_t) ((fBits & type) != 0);
   */
  Bool_t IsType(EProcessType type) { return TObject::TestBit(type); }

  /** @brief Function of initialization */
  virtual void Initialize() {}

  /**
   * @brief Function that compute length of interaction ?????
   * 
   * @param mat TGeoMaterial material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param lengths ?????
   * @param tid ?????
   */
  virtual void ComputeIntLen(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks,
                             Double_t *lengths, Int_t tid) = 0;

  /**
   * @brief Function that provides posterior steps
   * 
   * @param mat TGeoMaterial material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param nout ?????
   * @param tid ?????
   */
  virtual void PostStep(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t &nout,
                        Int_t tid) = 0;

  /**
   * @brief Post step type of intraction sampling function
   * @details Sampling:
   * 1. Target atom and type of the interaction for each primary tracks
   * 2. All inf. regarding sampling output is stored in the tracks
   * 
   * @param mat TGeoMaterial material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param tid  Track ID ?????
   */
  virtual void PostStepTypeOfIntrActSampling(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks,
                                             Int_t tid) = 0;

  /**
   * @brief Post step final state sampling function
   * @details Sampling final states for each primary tracks based on target atom and
   * interaction type sampled by PostStepTypeOfIntrActSampling;
   * updating primary track properties and inserting secondary tracks;
   * number of inserted secondary tracks will be stored in nout at termination;
   * 
   * @param mat TGeoMaterial material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param nout ?????
   * @param tid Track ID ?????
   */
  virtual void PostStepFinalStateSampling(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks,
                                          Int_t &nout, Int_t tid) = 0;

  /**
   * @todo  Need to be implemented
   */
  virtual void AtRest(Int_t /*ntracks*/, GeantTrack_v & /*tracks*/, Int_t & /*nout*/,
                      Int_t /*tid*/) {}
  GEANT_CUDA_DEVICE_CODE

  /**
   * @todo Need to be implemented
   */
  virtual void Eloss(TGeoMaterial * /*mat*/, Int_t /*ntracks*/, GeantTrack_v & /*tracks*/,
                     Int_t & /*nout*/, Int_t /*tid*/) {}
  GEANT_CUDA_DEVICE_CODE

  /**
   * @todo Need to be implemented
   */
  virtual void ApplyMsc(TGeoMaterial * /*mat*/, Int_t /*ntracks*/, GeantTrack_v & /*tracks*/,
                        Int_t /*tid*/) {}

  ClassDef(PhysicsProcess, 1) // Physics process base class
};

#endif
