//===--- ExN03ScoringData.h - Geant-V ----------------------- ---*- C++ -*-===//
//
//                     Geant-V Examples               
//
//===----------------------------------------------------------------------===//

#ifndef GEANT_ExN03_SCORING_DATA
#define GEANT_ExN03_SCORING_DATA

#include "GeantEvent.h"

//______________________________________________________________________________
//* @brief Class representing scoring data per layer for ExN03 (cumulated energy deposit and track length)
struct ExN03LayerDigit {
  float fEdepGap = 0.0f;   // Cumulated energy deposit in gap
  float fEdepAbs = 0.0f;   // Cumulated energy deposit in absorber
  float fLenGap = 0.0f;    // Cumulated step length in gap
  float fLenAbs = 0.0f;    // Cumulated step length in absorber
  
  ExN03LayerDigit &operator+=(const ExN03LayerDigit &other) {
    fEdepGap += other.fEdepGap;
    fEdepAbs += other.fEdepAbs;
    fLenGap += other.fLenGap;
    fLenAbs += other.fLenAbs;
    return *this;
  }
  
  inline
  void ScoreInGap(float edep, float len) { fEdepGap += edep; fLenGap += len; }

  inline
  void ScoreInAbs(float edep, float len) { fEdepAbs += edep; fLenAbs += len; }
  
  inline
  void Clear() { fEdepGap = fEdepAbs = fLenGap = fLenAbs = 0.0f; }
};
  

//______________________________________________________________________________
//* @brief Class representing scoring data per event for ExN03
class ExN03EventDigits {
private:
  int    fNlayers;     // Number of layers
  std::vector<ExN03LayerDigit> fDigits; //[fNlayers] Digits per event

public:
  ExN03EventDigits(int nlayers) : fNlayers(nlayers) {
    fDigits.resize(nlayers);
  }
  
  /** @brief Clear the data */
  void Clear() {
    for (auto digit : fDigits) digit.Clear();
  }
  
  /** @brief Getter for event digits */
  inline
  ExN03LayerDigit &GetDigit(int layer) { return fDigits[layer]; }
  
  /** @brief Merge with other digits */
  ExN03EventDigits &operator+=(const ExN03EventDigits &other) {
    for (int i=0; i<fNlayers; ++i)
      fDigits[i] += other.fDigits[i];
    return *this;
  }
  
  void Print() {
    for (int layer=0; layer<fNlayers; ++layer) {
      std::cout << "*** layer " << layer
                << ": edep_gap= " << fDigits[layer].fEdepGap << "   edep_abs= " << fDigits[layer].fEdepAbs
                << ": len_gap= " << fDigits[layer].fLenGap << "   len_abs= " << fDigits[layer].fLenAbs
                << std::endl;
    }
  }
};

//______________________________________________________________________________
//* @brief Class representing scoring data per event for ExN03
class ExN03ScoringData {
private:
  int fNbuff;  // Number of buffered events (should match GeantConfig::fNbuff
  int fNlayers;  // Number of layers
  std::vector<ExN03EventDigits> fDigits; // Digits data per event slot
  
public:
  ExN03ScoringData(int nbuffered, int nlayers) : fNbuff(nbuffered), fNlayers(nlayers) {
    fDigits.reserve(nbuffered);
    for (int i=0; i<nbuffered; ++i)
      fDigits.push_back(ExN03EventDigits(nlayers));
  }
  
  /** @brief Clear an event slot */
  void Clear(int evslot) {
    fDigits[evslot].Clear();
  }
  
  /** @brief Scoring data is summable */
  bool Merge(int evslot, const ExN03ScoringData &other) {
    fDigits[evslot] += other.GetDigits(evslot);
    return true;
  }
  /** @brief Print scoring data for event */
  void PrintDigits(Geant::GeantEvent *event) {
    fDigits[event->GetSlot()].Print();
  }

  /** @brief Getter for event digits */
  inline
  ExN03EventDigits &GetDigits(int evslot) { return fDigits[evslot]; }

  inline
  const ExN03EventDigits &GetDigits(int evslot) const { return fDigits[evslot]; }
};

#endif
