//===--- GeantOutput.h - Geant-V --------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantOutput.h
 * @brief Implementation of output of Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_OUTPUT
#define GEANT_OUTPUT

#include "globals.h"

/** @brief GeantOutput class */
class GeantOutput : public TObject {
public:
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

  double fCpuTime;       /** CPU time */
  Int_t fVolId;            /** Volume transporting this generation */
  Int_t fBasketGeneration; /** Current generation of baskets to be flushed */
  Int_t fGeneration;       /** Current generation for one basket */
  Int_t fNtracks;          /** Number of tracks in current generation */
  Int_t *fEvent;           //[fNtracks] /** Event */
  Int_t *fInd;             //[fNtracks] /** Track indices */
  Int_t *fProc;            //[fNtracks] /** Selected processes for each track */
  double *fX;            //[fNtracks] /** X positions */
  double *fY;            //[fNtracks] /** Y positions */
  double *fZ;            //[fNtracks] /** Z positions */
  double *fPx;           //[fNtracks] /** Px */
  double *fPy;           //[fNtracks] /** Py */
  double *fPz;           //[fNtracks] /** Pz */
  double *fE;            //[fNtracks] /** E */
  double *fPstep;        //[fNtracks] /** Physics step selected */
  double *fStep;         //[fNtracks] /** Current step */
  double *fSnext;        //[fNtracks] /** Snext distance */
  double *fSafety;       //[fNtracks] /** Snext distance */

private:

  /**
   * @brief GeantOutput copy constructor
   * @details Not allowed
   */
  GeantOutput(const GeantOutput &);
  
  /**
   * @brief Assignment operator
   * @details Not allowed
   */
  GeantOutput &operator=(const GeantOutput &); 

public:

  /** @brief GeantOutput constructor */
  GeantOutput()
      : TObject(), fCpuTime(0), fVolId(-1), fBasketGeneration(0), fGeneration(0), fNtracks(0),
        fEvent(0), fInd(0), fProc(0), fX(0), fY(0), fZ(0), fPx(0), fPy(0), fPz(0), fE(0), fPstep(0),
        fStep(0), fSnext(0), fSafety(0) {}

  /**
   * @brief GeantOutput destructor
   */
  virtual ~GeantOutput();

  /**
   * @brief Initialization function
   * 
   * @param size Size of output
   */
  void Init(Int_t size);

  /** @brief Reset function */
  void Reset();

  /**
   * @brief Set stamp function
   * 
   * @param volId Volume ID
   * @param basket_gen Number of basket od current generation
   * @param generation Nuber of generation
   * @param ntracks Number of tracks
   * @param cputime CPU time (by default 0)
   */
  void SetStamp(Int_t volId, Int_t basket_gen, Int_t generation, Int_t ntracks,
                double cputime = 0.) {
    fVolId = volId;
    fBasketGeneration = basket_gen;
    fGeneration = generation;
    fNtracks = ntracks;
    fCpuTime = cputime;
  }

  /**
   * @brief Set track function
   * 
   * @param ntrack Number of tracks
   * @param itrack Track
   * @param event Event
   * @param proc Selected processes for each track
   * @param x Z position
   * @param y Y position
   * @param z Z position
   * @param px Px
   * @param py Py
   * @param pz Pz
   * @param e Energy
   * @param pstep Physics step selected
   * @param step Current step
   * @param snext Next safety distance
   * @param safety Safety distance
   */
  void SetTrack(Int_t ntrack, Int_t itrack, Int_t event, Int_t proc, double x, double y,
                double z, double px, double py, double pz, double e, double pstep,
                double step, double snext, double safety);
  
  /**
   * @brief Function of setting tracks
   * 
   * @param ntrack Number of tracks
   * @param track Track to be setted
   */
  void SetTrack(Int_t ntrack, GeantTrack *track);

  ClassDef(GeantOutput, 1) // The transport output per generation
};
#endif
