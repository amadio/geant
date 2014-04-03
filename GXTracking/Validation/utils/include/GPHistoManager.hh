//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
//
// File: GPHistoManager.hh
//
// Purpose:
//   Central class to centralize management of histograms and similar objects.
//
// This is how to use it:
// 1. get a pointer:    GPHistoManager* hmgr = GPHistoManager::getInstance();
// 2. book a histogram: TH1F* hist = hmgr->book1F("name",100,0.1);
// 3. fill it (best):   hist->fill(value);
// 4. also works:       hmgr->fill("name",value);
//
// Last method can also be called from any other place, by accessing
// following steps 1 and 4.
//
// Note: make sure this class is destroyed at end of job, otherwise
// histograms may not be saved.
//
// 20131009 - Guilherme Lima - Created
//
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#ifndef _GPHISTOMANAGER_HH_
#define _GPHISTOMANAGER_HH_

#ifndef __CUDA_ARCH__
//#warning GPHistoManager.h: CUDA_ARCH undefined: Defining class GPHistoManager

// By making sure that GPHistoManager class is only valid if GPUPLOTS is defined,
// any attempt to use it under performance measurement conditions will cause compile time errors.
#ifdef GPUPLOTS

#include <string>
#include <map>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

// Map of SimCalorimeterHits keyed by cellID
typedef std::map<std::string,TH1F*> HistoMap;
typedef std::map<std::string,TH2F*> Histo2DMap;

class GPHistoManager {

public:
  static GPHistoManager& getInstance() {
    if( !_me ) _me = new GPHistoManager();
    return *_me;
  }

  void destroy();

  TH1F& getHisto(const std::string& name);

  // book new histograms and add to map, keyed by name
  TH1F& book1F(const std::string& name, int nbins, float xlo, float xhi) {
    TH1F* hist = new TH1F(name.c_str(),name.c_str(),nbins,xlo,xhi);
    _hists[name] = hist;
    return *_hists[name];
  }

  TH2F& book2F(const std::string& name, int nxbins, float xlo, float xhi, int nybins, float ylo, float yhi) {
    TH2F* hist = new TH2F(name.c_str(),name.c_str(),nxbins,xlo,xhi,nybins,ylo,yhi);
    _hists2D[name] = hist;
    return *_hists2D[name];
  }

  void fill(std::string& name, float value);
  void fill(std::string& name, float value, float val2);
  void fillFromVector(const std::string& name, std::vector<double>& values);
  void saveHistos();

  void bookBremHistos();
  void bookIoniHistos();
  void bookMscHistos();

private:
  // Constructor
  GPHistoManager();
  // Destructor
  ~GPHistoManager();

private:
  // ***** Member data  *****
  static GPHistoManager* _me;  // singleton pointer
  HistoMap _hists;            // histograms booked
  Histo2DMap _hists2D;          // histograms booked
  TFile _file;                // root file containing output histograms
};

#endif // GPUPLOTS

#else
//#warning GPHistoManager.h: Skipping  GPHistoManager definition due to CUDA_ARCH being defined
#endif // __CUDA_ARCH__

#endif // _GPHISTOMANAGER_HH_
