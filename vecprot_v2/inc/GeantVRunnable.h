//===--- GeantVRunnable.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantVRunnable.h
 * @brief Implementation of general interface to concurrent runnables. A runnable
 * is an object embedding data which is processed by the Run() method.
 * Running the object may produce a new runnable.
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_VRUNNABLE
#define GEANT_VRUNNABLE

/**
 * @brief Class GeantVRunnable
 * @details A general interface to concurrent runnables. A runnable
 * is an object embedding data which is processed by the Run()
 * method. Running the object may produce a new runnable.
 */
class GeantVRunnable {
public:
  
  /**
   * @brief GeantVRunnable constructor
   */
  GeantVRunnable();

  /**
   * @brief GeantVRunnable destructor
   */
  ~GeantVRunnable();

  /** @brief Run function */
  virtual GeantVRunnable *Run() = 0;

  ClassDef(GeantVRunnable, 1) // ABC for concurrent runnables
};
