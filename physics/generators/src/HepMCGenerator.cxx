#include "HepMCGenerator.h"
#include "base/Global.h"
using vecgeom::kPi;
#include "Geant/Typedefs.h"
#include "Geant/Error.h"

#include "TaskData.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/ReaderAscii.h"
#ifdef USE_ROOT
#include "HepMC/ReaderRoot.h"
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
HepMCGenerator::HepMCGenerator()
  : input_file(0), search(0) {}

//______________________________________________________________________________
HepMCGenerator::HepMCGenerator(std::string &filename) : input_file(0), search(0) {
  if (filename.substr(filename.find_last_of(".") + 1) == "hepmc3") {
    input_file = new HepMC::ReaderAscii(filename);
  }
#ifdef USE_ROOT
  else if (filename.substr(filename.find_last_of(".") + 1) == "root") {
    input_file = new HepMC::ReaderRoot(filename);
  }
#endif
  else {
    std::cout << "Unrecognized filename extension (must be .hepmc3 or .root)" << std::endl;
  }
  Particle_t::CreateParticles();
}

//______________________________________________________________________________
HepMCGenerator::~HepMCGenerator() {
  delete input_file;
  delete search;
}

//______________________________________________________________________________
void HepMCGenerator::InitPrimaryGenerator() {}

//______________________________________________________________________________
EventInfo HepMCGenerator::NextEvent(geant::TaskData* /*td*/) {
  //
  // Delete previous event
  delete search;

  HepMC::GenEvent evt(HepMC::Units::GEV, HepMC::Units::MM);

  if (!(input_file->read_event(evt)))
    geant::Fatal("HepMCGenerator::ImportTracks", "No more particles to read!");

  //  std::cout << std::endl
  //            << "Find all stable particles: " << std::endl;

  search = new HepMC::FindParticles(evt, HepMC::FIND_ALL, HepMC::STATUS == 1); // && HepMC::STATUS_SUBCODE == 0);

  int ntracks = 0;
  double eta, phi, pmom = 0;
  for (const HepMC::GenParticlePtr &genpart : search->results()) {
    if (fEtaCut || fMomCut)
      pmom = sqrt(genpart->momentum().px() * genpart->momentum().px() +
                  genpart->momentum().py() * genpart->momentum().py() +
                  genpart->momentum().pz() * genpart->momentum().pz());
    // genpart->print();
    if (fEtaCut) {
      if (pmom == genpart->momentum().pz())
        eta = 1.E30;
      else
        eta = 0.5 * log((pmom + genpart->momentum().pz()) / (pmom - genpart->momentum().pz()));
      if (eta < fEtaMin || eta > fEtaMax)
        continue;
    }
    if (fPhiCut) {
      // Phi in 0,2pi
      phi = kPi + atan2(-genpart->momentum().py(), -genpart->momentum().px());
      if (phi < fPhiMin || phi > fPhiMax)
        continue;
    }
    if (fMomCut) {
      if (pmom < fPMin || pmom > fPMax)
        continue;
    }
    ntracks++;
  }

  int ntot = search->results().size();

  std::cout << std::endl
            << "Number of stable particles ";
  if (fEtaCut || fPhiCut)
    std::cout << "within selected cuts: ";
  std::cout << ntracks;
  if (ntot > ntracks)
    std::cout << " out of " << ntot;
  std::cout << std::endl;

  EventInfo current;
  current.ntracks = ntracks;
  current.xvert = evt.event_pos().x();
  current.yvert = evt.event_pos().y();
  current.zvert = evt.event_pos().z();
  current.tvert = evt.event_pos().t();
  return current;
}

//______________________________________________________________________________
void HepMCGenerator::GetTrack(int n, geant::Track &gtrack, geant::TaskData* /*td*/) {
  //  const HepMC::GenParticlePtr &genpart = search->results()[n];
  int itr = 0;
  double eta, phi, pmom = 0;
  for (const HepMC::GenParticlePtr &genpart : search->results()) {
    if (itr++ < n)
      continue;
    if (fEtaCut || fMomCut || fPhiCut) {
      pmom = sqrt(genpart->momentum().px() * genpart->momentum().px() +
                  genpart->momentum().py() * genpart->momentum().py() +
                  genpart->momentum().pz() * genpart->momentum().pz());
      if (fEtaCut) {
        if (pmom == genpart->momentum().pz())
          eta = 1.E30;
        else
          eta = 0.5 * log((pmom + genpart->momentum().pz()) / (pmom - genpart->momentum().pz()));
        if (eta < fEtaMin || eta > fEtaMax)
          continue;
      }
      if (fPhiCut) {
        // Phi in 0,2pi
        phi = kPi + atan2(-genpart->momentum().py(), -genpart->momentum().px());
        if (phi < fPhiMin || phi > fPhiMax)
          continue;
      }
      if (fMomCut) {
        if (pmom < fPMin || pmom > fPMax)
          continue;
      }
    }

    // here I have to create Tracks
    int pdg = genpart->pid();
    gtrack.SetPDG(pdg);
//    gtrack.SetPDG(2212);

    gtrack.SetGVcode(TPartIndex::I()->PartIndex(gtrack.PDG()));
    const Particle_t *const &part = &Particle_t::GetParticle(gtrack.PDG());
    gtrack.SetCharge(part->Charge());
    gtrack.SetMass(part->Mass());

    if ((bool)genpart->production_vertex()) {
      gtrack.SetPosition(genpart->production_vertex()->position().x(),
                         genpart->production_vertex()->position().y(),
                         genpart->production_vertex()->position().z());
    } else {
      gtrack.SetPosition(0., 0., 0.);
    }


    // Compute momentum from energy/mass
    //double p = sqrt((gtrack.E() - gtrack.Mass()) * (gtrack.E() + gtrack.Mass()));
    // Momentum from generator
    double ptrack =
        sqrt(genpart->momentum().px() * genpart->momentum().px() + genpart->momentum().py() * genpart->momentum().py() +
             genpart->momentum().pz() * genpart->momentum().pz());

    gtrack.SetP(ptrack);
    gtrack.SetE(sqrt(ptrack*ptrack+gtrack.Mass()*gtrack.Mass()));
//    std::cout << "track momentum = " << ptrack*1000. << std::endl;
    // Correctly normalize direction
    gtrack.SetDirection(genpart->momentum().px() / ptrack,
                        genpart->momentum().py() / ptrack,
                        genpart->momentum().pz() / ptrack);
    return;
  }
}

//______________________________________________________________________________
void HepMCGenerator::GetTrack(int n, double &tpx, double &tpy, double &tpz, double &te, double &x0, double &y0,
                              double &z0, int &pdg) {

  const HepMC::GenParticlePtr &genpart = search->results()[n];
  // here I have to create Tracks
  pdg = genpart->pid();
//  pdg = 2212;
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
    const Particle_t *const &part = &Particle_t::GetParticle(pdg);
  te = sqrt(tpx*tpx+tpy*tpy+tpz*tpz+part->Mass()*part->Mass()*1E6);
//  std::cout << "track momentum = " << te << std::endl;
//  te = genpart->momentum().e() * 1000.0;

  /*
  //  const HepMC::GenParticlePtr &genpart = search->results()[n];
    int itr = 0;
    double eta, phi, pmom=0;
    for (const HepMC::GenParticlePtr &genpart : search->results()) {
      if (fEtaCut || fMomCut)
        pmom = sqrt(genpart->momentum().px()*genpart->momentum().px() +
                           genpart->momentum().py()*genpart->momentum().py() +
                           genpart->momentum().pz()*genpart->momentum().pz());
      if (fEtaCut) {
        pmom = sqrt(genpart->momentum().px()*genpart->momentum().px() +
                                 genpart->momentum().py()*genpart->momentum().py() +
                                 genpart->momentum().pz()*genpart->momentum().pz());
        if (pmom == genpart->momentum().pz()) eta = 1.E30;
        else eta = 0.5*log((pmom+genpart->momentum().pz())/(pmom-genpart->momentum().pz()));
        if (eta<fEtaMin || eta>fEtaMax) continue;
      }
      if (fPhiCut) {
        // Phi in 0,2pi
        phi = kPi+atan2(-genpart->momentum().py(),-genpart->momentum().px());
        if (phi<fPhiMin || phi>fPhiMax) continue;
      }
      if (fMomCut) {
        if (pmom<fPMin || pmom>fPMax) continue;
      }
      if (itr++<n) continue;
      // here I have to create Tracks
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

} // GEANT_IMPL_NAMESPACE
} // Geant
