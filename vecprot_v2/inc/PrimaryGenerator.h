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
#include "TNamed.h"

#include "GeantFwd.h"

using vecgeom::kRadToDeg;
using vecgeom::kDegToRad;

/**
 * @brief Class of primary generators
 */
class PrimaryGenerator : public TNamed {
protected:
  Bool_t fEtaCut;   // flag for specifying a cut in eta
  Bool_t fPhiCut;   // flag for specifying a cut in phi
  Bool_t fMomCut;   // flag for specifying a cut in momentum
  Double_t fEtaMin; // minimum eta
  Double_t fEtaMax; // maximum eta
  Double_t fPhiMin; // minimum phi
  Double_t fPhiMax; // maximum phi
  Double_t fPMin;   // minimum momentum
  Double_t fPMax;   // maximum momentum
public:
  PrimaryGenerator()
      : TNamed(), fEtaCut(false), fPhiCut(false), fMomCut(false), fEtaMin(0), fEtaMax(0), fPhiMin(0), fPhiMax(0),
        fPMin(0), fPMax(0) {}
  virtual ~PrimaryGenerator() {}
  static Double_t EtaToTheta(Double_t eta) { return (2. * atan(exp(-eta)) * kRadToDeg); }
  static Double_t ThetaToEta(Double_t theta) { return (-log(tan(0.5 * theta * kDegToRad))); }

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
   * @param n Track index
   * @param gtrack track
   */
  virtual void GetTrack(Int_t n, Geant::GeantTrack &gtrack) = 0;

  /** @brief Getter for eta cut flag */
  Bool_t HasEtaCut() const { return fEtaCut; }

  /** @brief Getter for phi cut flag */
  Bool_t HasPhiCut() const { return fPhiCut; }

  /** @brief Getter for momentum cut flag */
  Bool_t HasMomCut() const { return fMomCut; }

  /** @brief Setter for user eta range */
  void SetEtaRange(Double_t etamin, Double_t etamax) {
    fEtaCut = true;
    fEtaMin = etamin;
    fEtaMax = etamax;
  }

  /** @brief Getter for user eta range */
  void GetEtaRange(Double_t &etamin, Double_t &etamax) const {
    etamin = fEtaMin;
    etamax = fEtaMax;
  }

  /** @brief Setter for user theta range */
  void SetThetaRange(Double_t thetamin_deg, Double_t thetamax_deg) {
    fEtaCut = true;
    fEtaMin = ThetaToEta(thetamax_deg);
    fEtaMax = ThetaToEta(thetamin_deg);
  }

  /** @brief Getter for user theta range */
  void GetThetaRange(Double_t &thetamin_deg, Double_t &thetamax_deg) const {
    thetamin_deg = EtaToTheta(fEtaMax);
    thetamax_deg = EtaToTheta(fEtaMin);
  }

  /** @brief Setter for user phi range */
  void SetPhiRange(Double_t phimin_deg, Double_t phimax_deg) {
    fPhiCut = true;
    fPhiMin = phimin_deg;
    fPhiMax = phimax_deg;
  }

  /** @brief Getter for user phi range */
  void GetPhiRange(Double_t &phimin, Double_t &phimax) const {
    phimin = fPhiMin;
    phimax = fPhiMax;
  }

  /** @brief Setter for user momentum range */
  void SetMomRange(Double_t pmin, Double_t pmax) {
    fMomCut = true;
    fPMin = pmin;
    fPMax = pmax;
  }

  /** @brief Getter for user momentum range */
  void GetMomRange(Double_t &pmin, Double_t &pmax) const {
    pmin = fPMin;
    pmax = fPMax;
  }
};

#endif
