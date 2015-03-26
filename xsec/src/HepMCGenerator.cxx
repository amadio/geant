#include "HepMCGenerator.h"

#include "TMath.h"
#include "TDatabasePDG.h"
#include "GeantTrack.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IO/IO_GenEvent.h"
#include "HepMC/IO/IO_Root.h"

ClassImp(HepMCGenerator)

    //______________________________________________________________________________
    HepMCGenerator::HepMCGenerator()
    : input_file(0), search(0) {}

//______________________________________________________________________________
HepMCGenerator::HepMCGenerator(std::string &filename) : input_file(0), search(0) {
  if (filename.substr(filename.find_last_of(".") + 1) == "hepmc3") {
    input_file = new HepMC::IO_GenEvent(filename, std::ios::in);
  } else if (filename.substr(filename.find_last_of(".") + 1) == "root") {
    input_file = new HepMC::IO_Root(filename, std::ios::in);
  } else {
    std::cout << "Unrecognized filename extension (must be .hepmc3 or .root)" << std::endl;
  }
}

//______________________________________________________________________________
HepMCGenerator::~HepMCGenerator() {
  delete input_file;
  delete search;
}

//______________________________________________________________________________
void HepMCGenerator::InitPrimaryGenerator() {}

//______________________________________________________________________________
Int_t HepMCGenerator::NextEvent() {
  //
  // Delete previous event
  delete search;

  HepMC::GenEvent evt(HepMC::Units::GEV, HepMC::Units::MM);

  if (!(input_file->fill_next_event(evt)))
    Fatal("ImportTracks", "No more particles to read!");

//  std::cout << std::endl
//            << "Find all stable particles: " << std::endl;

  search = new HepMC::FindParticles(evt, HepMC::FIND_ALL, HepMC::STATUS == 1 && HepMC::STATUS_SUBCODE == 0);

  for (const HepMC::GenParticlePtr &genpart : search->results()) {
//    genpart->print();
  }

  Int_t ntracks = search->results().size();

  std::cout << std::endl
            << "Number of stable particles " << ntracks << std::endl;
//ntracks = 1;
  return ntracks;
}

//______________________________________________________________________________
void HepMCGenerator::GetTrack(Int_t n, GeantTrack &gtrack) {
  const HepMC::GenParticlePtr &genpart = search->results()[n];

  // here I have to create GeantTracks
  int pdg = genpart->pdg_id();

  gtrack.SetPDG(pdg);

  gtrack.SetG5code(TPartIndex::I()->PartIndex(pdg));
  TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(gtrack.fPDG);

  gtrack.SetCharge(part->Charge() / 3.);
  gtrack.SetMass(part->Mass());

  if ((bool)genpart->production_vertex()) {
    gtrack.fXpos = genpart->production_vertex()->position().x();
    gtrack.fYpos = genpart->production_vertex()->position().y();
    gtrack.fZpos = genpart->production_vertex()->position().z();
  } else {
    gtrack.fXpos = 0;
    gtrack.fYpos = 0;
    gtrack.fZpos = 0;
  }

  gtrack.fE = genpart->momentum().e(); // e- 30MeV

  Double_t p = TMath::Sqrt((gtrack.E() - gtrack.Mass()) * (gtrack.E() + gtrack.Mass()));

  gtrack.SetP(p);
  gtrack.fXdir = genpart->momentum().px() / p;
  gtrack.fYdir = genpart->momentum().py() / p;
  gtrack.fZdir = genpart->momentum().pz() / p;

}

//______________________________________________________________________________
void HepMCGenerator::GetTrack(Int_t n, Double_t &tpx, Double_t &tpy, Double_t &tpz, Double_t &te,
                              Double_t &x0, Double_t &y0, Double_t &z0, Int_t &pdg) {
  const HepMC::GenParticlePtr &genpart = search->results()[n];
  // here I have to create GeantTracks
  pdg = genpart->pdg_id();
  if ((bool)genpart->production_vertex()) {
    // current default unit is [mm] that is the default Geant4 length unit as well
    x0 = genpart->production_vertex()->position().x();
    y0 = genpart->production_vertex()->position().y();
    z0 = genpart->production_vertex()->position().z();
  } else {
    x0 = 0.0;
    y0 = 0.0;
    z0 = 0.0;
  }

  // current default unit is [GeV] so change it to [MeV] for default Geant4 energy unit 
  tpx = genpart->momentum().px()*1000.0;
  tpy = genpart->momentum().py()*1000.0;
  tpz = genpart->momentum().pz()*1000.0;     
  te  = genpart->momentum().e() *1000.0;

}

