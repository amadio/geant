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

#include "base/Global.h"

#include "GeantFwd.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

using vecgeom::kRadToDeg;
using vecgeom::kDegToRad;

class GeantTaskData;

struct EventInfo {
  int      ntracks;    // number of tracks
  double   xvert;      // x position
  double   yvert;      // y position
  double   zvert;      // z position
  double   tvert;      // time
// ... to be extended
  EventInfo() : xvert(0), yvert(0), zvert(0) {}
};

/**
 * @brief Class of primary generators
 */
class PrimaryGenerator {
protected:
  bool fEtaCut;   // flag for specifying a cut in eta
  bool fPhiCut;   // flag for specifying a cut in phi
  bool fMomCut;   // flag for specifying a cut in momentum
  double fEtaMin; // minimum eta
  double fEtaMax; // maximum eta
  double fPhiMin; // minimum phi
  double fPhiMax; // maximum phi
  double fPMin;   // minimum momentum
  double fPMax;   // maximum momentum
public:
  PrimaryGenerator()
      : fEtaCut(false), fPhiCut(false), fMomCut(false), fEtaMin(0), fEtaMax(0), fPhiMin(0), fPhiMax(0),
        fPMin(0), fPMax(0) {}
  virtual ~PrimaryGenerator() {}
  static double EtaToTheta(double eta) { return (2. * atan(exp(-eta)) * kRadToDeg); }
  static double ThetaToEta(double theta) { return (-log(tan(0.5 * theta * kDegToRad))); }

  /**
   * @brief Pure virtual function of initialization of primary generator
   * @details Set one Track primary track properties
   */
  virtual void InitPrimaryGenerator() = 0;

  /** @brief  Pure virtual function that produce next event
    *
    * @param td thread local data pointer
    */
  virtual EventInfo NextEvent(geant::GeantTaskData* td) = 0;

  /**
   * @brief Pure virtual function that returns track
   *
   * @param n Track index
   * @param gtrack track
   * @param td thread local data pointer
   */
  virtual void GetTrack(int n, geant::Track &gtrack, geant::GeantTaskData* td) = 0;

  /** @brief Getter for eta cut flag */
  bool HasEtaCut() const { return fEtaCut; }

  /** @brief Getter for phi cut flag */
  bool HasPhiCut() const { return fPhiCut; }

  /** @brief Getter for momentum cut flag */
  bool HasMomCut() const { return fMomCut; }

  /** @brief Setter for user eta range */
  void SetEtaRange(double etamin, double etamax) {
    fEtaCut = true;
    fEtaMin = etamin;
    fEtaMax = etamax;
  }

  /** @brief Getter for user eta range */
  void GetEtaRange(double &etamin, double &etamax) const {
    etamin = fEtaMin;
    etamax = fEtaMax;
  }

  /** @brief Setter for user theta range */
  void SetThetaRange(double thetamin_deg, double thetamax_deg) {
    fEtaCut = true;
    fEtaMin = ThetaToEta(thetamax_deg);
    fEtaMax = ThetaToEta(thetamin_deg);
  }

  /** @brief Getter for user theta range */
  void GetThetaRange(double &thetamin_deg, double &thetamax_deg) const {
    thetamin_deg = EtaToTheta(fEtaMax);
    thetamax_deg = EtaToTheta(fEtaMin);
  }

  /** @brief Setter for user phi range */
  void SetPhiRange(double phimin_deg, double phimax_deg) {
    fPhiCut = true;
    fPhiMin = phimin_deg;
    fPhiMax = phimax_deg;
  }

  /** @brief Getter for user phi range */
  void GetPhiRange(double &phimin, double &phimax) const {
    phimin = fPhiMin;
    phimax = fPhiMax;
  }

  /** @brief Setter for user momentum range */
  void SetMomRange(double pmin, double pmax) {
    fMomCut = true;
    fPMin = pmin;
    fPMax = pmax;
  }

  /** @brief Getter for user momentum range */
  void GetMomRange(double &pmin, double &pmax) const {
    pmin = fPMin;
    pmax = fPMax;
  }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
