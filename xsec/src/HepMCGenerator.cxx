#include "HepMCGenerator.h"

#include "TMath.h"
#include "TDatabasePDG.h"
#include "GeantTrack.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/ReaderAscii.h"
#include "HepMC/ReaderRoot.h"

ClassImp(HepMCGenerator)

    //______________________________________________________________________________
    HepMCGenerator::HepMCGenerator()
    : input_file(0), search(0) {}

//______________________________________________________________________________
HepMCGenerator::HepMCGenerator(std::string &filename) : input_file(0), search(0) {
  if (filename.substr(filename.find_last_of(".") + 1) == "hepmc3") {
    input_file = new HepMC::ReaderAscii(filename);
  } else if (filename.substr(filename.find_last_of(".") + 1) == "root") {
    input_file = new HepMC::ReaderRoot(filename);
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

  if (!(input_file->read_event(evt)))
    Fatal("ImportTracks", "No more particles to read!");

  //  std::cout << std::endl
  //            << "Find all stable particles: " << std::endl;

  search = new HepMC::FindParticles(evt, HepMC::FIND_ALL, HepMC::STATUS == 1); // && HepMC::STATUS_SUBCODE == 0);

  Int_t ntracks = 0;
  Double_t eta, phi, pmom = 0;
  for (const HepMC::GenParticlePtr &genpart : search->results()) {
    if (fEtaCut || fMomCut)
      pmom = TMath::Sqrt(genpart->momentum().px() * genpart->momentum().px() +
                         genpart->momentum().py() * genpart->momentum().py() +
                         genpart->momentum().pz() * genpart->momentum().pz());
    // genpart->print();
    if (fEtaCut) {
      if (pmom == genpart->momentum().pz())
        eta = 1.E30;
      else
        eta = 0.5 * TMath::Log((pmom + genpart->momentum().pz()) / (pmom - genpart->momentum().pz()));
      if (eta < fEtaMin || eta > fEtaMax)
        continue;
    }
    if (fPhiCut) {
      // Phi in 0,2pi
      phi = TMath::Pi() + TMath::ATan2(-genpart->momentum().py(), -genpart->momentum().px());
      if (phi < fPhiMin || phi > fPhiMax)
        continue;
    }
    if (fMomCut) {
      if (pmom < fPMin || pmom > fPMax)
        continue;
    }
    ntracks++;
  }

  Int_t ntot = search->results().size();

  std::cout << std::endl
            << "Number of stable particles ";
  if (fEtaCut || fPhiCut)
    std::cout << "within selected cuts: ";
  std::cout << ntracks;
  if (ntot > ntracks)
    std::cout << " out of " << ntot;
  std::cout << std::endl;

  return ntracks;
}

//______________________________________________________________________________
void HepMCGenerator::GetTrack(Int_t n, Geant::GeantTrack &gtrack) {
  //  const HepMC::GenParticlePtr &genpart = search->results()[n];
  Int_t itr = 0;
  Double_t eta, phi, pmom = 0;
  for (const HepMC::GenParticlePtr &genpart : search->results()) {
    if (fEtaCut || fMomCut)
      pmom = TMath::Sqrt(genpart->momentum().px() * genpart->momentum().px() +
                         genpart->momentum().py() * genpart->momentum().py() +
                         genpart->momentum().pz() * genpart->momentum().pz());
    if (fEtaCut) {
      pmom = TMath::Sqrt(genpart->momentum().px() * genpart->momentum().px() +
                         genpart->momentum().py() * genpart->momentum().py() +
                         genpart->momentum().pz() * genpart->momentum().pz());
      if (pmom == genpart->momentum().pz())
        eta = 1.E30;
      else
        eta = 0.5 * TMath::Log((pmom + genpart->momentum().pz()) / (pmom - genpart->momentum().pz()));
      if (eta < fEtaMin || eta > fEtaMax)
        continue;
    }
    if (fPhiCut) {
      // Phi in 0,2pi
      phi = TMath::Pi() + TMath::ATan2(-genpart->momentum().py(), -genpart->momentum().px());
      if (phi < fPhiMin || phi > fPhiMax)
        continue;
    }
    if (fMomCut) {
      if (pmom < fPMin || pmom > fPMax)
        continue;
    }
    if (itr++ < n)
      continue;

    // here I have to create GeantTracks
    int pdg = genpart->pid();
    gtrack.SetPDG(pdg);

    gtrack.SetGVcode(TPartIndex::I()->PartIndex(pdg));
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

    // Compute momentum from energy/mass
    Double_t p = TMath::Sqrt((gtrack.E() - gtrack.Mass()) * (gtrack.E() + gtrack.Mass()));
    // Momentum from generator
    Double_t ptrack = TMath::Sqrt(genpart->momentum().px() * genpart->momentum().px() +
                                  genpart->momentum().py() * genpart->momentum().py() +
                                  genpart->momentum().pz() * genpart->momentum().pz());

    gtrack.SetP(p);
    // Correctly normalize direction
    gtrack.fXdir = genpart->momentum().px() / ptrack;
    gtrack.fYdir = genpart->momentum().py() / ptrack;
    gtrack.fZdir = genpart->momentum().pz() / ptrack;
    return;
  }
}

//______________________________________________________________________________
void HepMCGenerator::GetTrack(Int_t n, Double_t &tpx, Double_t &tpy, Double_t &tpz, Double_t &te, Double_t &x0,
                              Double_t &y0, Double_t &z0, Int_t &pdg) {

  const HepMC::GenParticlePtr &genpart = search->results()[n];
  // here I have to create GeantTracks
  pdg = genpart->pid();
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
  tpx = genpart->momentum().px() * 1000.0;
  tpy = genpart->momentum().py() * 1000.0;
  tpz = genpart->momentum().pz() * 1000.0;
  te = genpart->momentum().e() * 1000.0;

  /*
  //  const HepMC::GenParticlePtr &genpart = search->results()[n];
    Int_t itr = 0;
    Double_t eta, phi, pmom=0;
    for (const HepMC::GenParticlePtr &genpart : search->results()) {
      if (fEtaCut || fMomCut)
        pmom = TMath::Sqrt(genpart->momentum().px()*genpart->momentum().px() +
                           genpart->momentum().py()*genpart->momentum().py() +
                           genpart->momentum().pz()*genpart->momentum().pz());
      if (fEtaCut) {
        pmom = TMath::Sqrt(genpart->momentum().px()*genpart->momentum().px() +
                                 genpart->momentum().py()*genpart->momentum().py() +
                                 genpart->momentum().pz()*genpart->momentum().pz());
        if (pmom == genpart->momentum().pz()) eta = 1.E30;
        else eta = 0.5*TMath::Log((pmom+genpart->momentum().pz())/(pmom-genpart->momentum().pz()));
        if (eta<fEtaMin || eta>fEtaMax) continue;
      }
      if (fPhiCut) {
        // Phi in 0,2pi
        phi = TMath::Pi()+TMath::ATan2(-genpart->momentum().py(),-genpart->momentum().px());
        if (phi<fPhiMin || phi>fPhiMax) continue;
      }
      if (fMomCut) {
        if (pmom<fPMin || pmom>fPMax) continue;
      }
      if (itr++<n) continue;
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
  */
}
