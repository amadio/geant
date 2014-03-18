
#ifndef _GPHISTOMANAGER_CC_
#define _GPHISTOMANAGER_CC_
//#warning GPHistoManager.cc: implementing GPHistoManager methods
#include "include/GPHistoManager.hh"
#include <iostream>

using std::string;
using std::vector;
using std::pair;
using std::cout;
using std::endl;

GPHistoManager* GPHistoManager::_me = NULL;

// To be used instead of constructor
//static GPHistoManager& GPHistoManager::getInstance() {


GPHistoManager::GPHistoManager() : _file("histos.root","recreate") { }
GPHistoManager::~GPHistoManager() {
  printf("Destructor: Writing (header)...\n");
  _file.Write();
  printf("Destructor: Closing...\n");
  _file.Close();
  printf("Destructor: Done!\n");
}


void GPHistoManager::destroy() {
  // first save and delete all histograms
  for(HistoMap::iterator it=_hists.begin(); it!=_hists.end(); ++it) {
    TH1F* hist = it->second;
    hist->Write();
    delete hist;
    _hists.erase(it);
  }
  // then destroy itself
  if(_me) delete _me;
}

//TH1F& GPHistoManager::book1F(const string& name, int nbins, float xlo, float xhi) {

TH1F& GPHistoManager::getHisto(const std::string& name) {
  // if it does not exist, create it now
  TH1F* phist = _hists[name];
  if(!phist) this->book1F(name,50,0,100);  // default, if not explicitly booked
  return *_hists[name];
}

void GPHistoManager::fill(string& name, float value) {
  _hists[name]->Fill(value);
}

void GPHistoManager::fill(string& name, float value, float val2) {
  _hists2D[name]->Fill(value,val2);
}

void GPHistoManager::fillFromVector(const std::string& name, vector<double>& values) {
  TH1F& hist = getHisto(name);
  vector<double>::const_iterator iter;
  for(iter=values.begin(); iter!=values.end(); ++iter) {
    hist.Fill( *iter );
  }
  // clear vector after values have been used
  values.clear();
}

void GPHistoManager::saveHistos() {
  for(HistoMap::iterator it=_hists.begin(); it!=_hists.end(); ++it) {
    it->second->Write();
  }

  for(Histo2DMap::iterator it=_hists2D.begin(); it!=_hists2D.end(); ++it) {
    it->second->Write();
  }
}

void GPHistoManager::bookBremHistos() {

  this->book1F("hbremAngle",100,0,0.2).Sumw2();
  this->book1F("hbremEnergy",100,0,1000).Sumw2();
  this->book1F("hbremSecAngle",100,0,0.2).Sumw2();
  this->book1F("hbremSecEnergy",103,-10,1020).Sumw2();

  this->book1F("hbremEnergyLoss",100,0,10).Sumw2();
  this->book1F("hbremStepLengthPost",100,0,10).Sumw2();
  this->book1F("hbremStepLengthAlong",100,0,10).Sumw2();
  this->book1F("hbremStepLengthDiff",100,0,10).Sumw2();

  this->book1F("hbremPreStepLambda",100,0,2).Sumw2();
  this->book1F("hbremNbOfIntLengthLeft",110,-1,10).Sumw2();
  this->book1F("hbremPreStepScaledEnergy",103,-10,1020).Sumw2();
}

void GPHistoManager::bookIoniHistos() {

  this->book1F("hioniAngle",100,0,0.5).Sumw2();
  this->book1F("hioniEnergy",102,-10,1020).Sumw2();
  this->book1F("hioniSecAngle",100,0,1).Sumw2();
  this->book1F("hioniSecEnergy",100,0,500).Sumw2();

  this->book1F("hioniEnergyLoss",100,0,100).Sumw2();
  this->book1F("hioniStepLengthPost",100,0,300).Sumw2();
  this->book1F("hioniStepLengthAlong",100,0,300).Sumw2();
  this->book1F("hioniStepLengthDiff",100,0,50).Sumw2();

  this->book1F("hioniPreStepLambda",100,0.04,0.05).Sumw2();
  this->book1F("hioniNbOfIntLengthLeft",110,-1,15).Sumw2();
  this->book1F("hioniPreStepScaledEnergy",103,-10,1020).Sumw2();
  this->book1F("hioniMassRatio",50,0.998,1.002).Sumw2();

  this->book1F("hioniDedxForScaledEnergyTimesLength",100,0,20).Sumw2();
  this->book1F("hioniElossFromKinEnergyMinusScaledEnergyForLoss",100,0,20).Sumw2();
  this->book1F("hioniElossFromSampleFluctuations",100,0,20).Sumw2();
  this->book1F("hioniEloss",100,0,20).Sumw2();
}

void GPHistoManager::bookMscHistos() {

  this->book1F("hmscAngle",100,0,10);
  this->book1F("hmscEnergy",102,-10,1010);

  this->book1F("hmscEnergyLoss",100,0,10);
  this->book1F("hmscStepLength",100,0,10);

  this->book1F("hmscPreStepLambda",110,-1,10);
  this->book1F("hmscNbOfIntLengthLeft",110,-1,10);
  this->book1F("hmscPreStepScaledEnergy",101,-10,1000);
}

#else
#warning GPHistoManager.cc: Skipping member function implementations
#endif // _GPHISTOMANAGER_CC_
