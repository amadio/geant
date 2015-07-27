//===--- ExN03Application.h - Geant-V ------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file ExN03Application.h
 * @brief Implementation of Geant4 example N03 application for Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_ExN03Application
#define GEANT_ExN03Application

#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif

#include "Geant/Typedefs.h"

#include "GeantFwd.h"

/** @brief ExN03Application class */
class ExN03Application : public GeantVApplication {
  static const int kNlayers = 10;
  static const int kMaxThreads = 36;

private:
  Bool_t fInitialized;                       /** Initialized flag */
  int fIdGap;                              /** ID for the gap volume */
  int fIdAbs;                              /** ID for the absorber volume */
  Float_t fEdepGap[kNlayers][kMaxThreads];   /** Energy deposition per layer */
  Float_t fLengthGap[kNlayers][kMaxThreads]; /** step length in every layer */
  Float_t fEdepAbs[kNlayers][kMaxThreads];   /** Energy deposition per layer */
  Float_t fLengthAbs[kNlayers][kMaxThreads]; /** Step length in every layer */
  GeantFactory<MyHit> *fFactory;             /** Hits factory */
  
  /**
   * @brief Copy constructor ExN03Application
   * * @todo Still not implemented
   */
  ExN03Application(const ExN03Application &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  ExN03Application &operator=(const ExN03Application &);
public:

  /** @brief Constructor ExN03Application */
  ExN03Application();

  /** @brief Destructor ExN03Application */
  virtual ~ExN03Application() {}

  /**
   * @brief Function of initialization
   */
  virtual Bool_t Initialize();

  /**
   * @brief Function that provides step manager 
   * 
   * @param tid ?????
   * @param npart ?????
   * @param tracks GeantV tracks
   */
  virtual void StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td);

  /**
   * @brief Function of digitization
   * 
   * @param event Event that should be digitized
   */
  virtual void Digitize(int event);

  /** @brief User FinishRun function */
  virtual void FinishRun() {}

  ClassDef(ExN03Application, 1) // User application
};
#endif
