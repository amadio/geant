#include "HepMCGenerator.h"

#include "TMath.h"
#include "TDatabasePDG.h"
#include "GeantTrack.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"


ClassImp(HepMCGenerator)

//______________________________________________________________________________
HepMCGenerator::HepMCGenerator(): input_file(0)
{
}


//______________________________________________________________________________
HepMCGenerator::HepMCGenerator(std::string& filename)
{
    input_file = new HepMC::IO_GenEvent(filename, std::ios::in);
}

//______________________________________________________________________________
HepMCGenerator::~HepMCGenerator()
{
    input_file->close();

    delete input_file;
}


//______________________________________________________________________________
void HepMCGenerator::InitPrimaryGenerator(){
}

//______________________________________________________________________________
Int_t HepMCGenerator::NextEvent(){
    //
    HepMC::GenEvent evt(HepMC::Units::GEV,HepMC::Units::MM);
    
    input_file->fill_next_event(evt);
    
    if(input_file->rdstate()) Fatal("ImportTracks","No more particles to read!");
    
    std::cout<<std::endl<<"Find all stable particles: "<<std::endl;

    search = new HepMC::FindParticles(evt, HepMC::FIND_ALL, HepMC::STATUS == 1 && HepMC::STATUS_SUBCODE == 0);
    
    for( const HepMC::GenParticlePtr &genpart: search->results() ) {
      genpart->print();
    }  

    Int_t ntracks =  search->results().size();

    std::cout << std::endl << "Number of stable particles " << ntracks << std::endl;

    return ntracks;
}

//______________________________________________________________________________
void HepMCGenerator::GetTrack(Int_t n, GeantTrack &gtrack){

  const HepMC::GenParticlePtr &genpart = search->results()[n]; 

  // here I have to create GeantTracks
  int pdg = genpart->pdg_id();

  gtrack.SetPDG(pdg);

  gtrack.SetG5code(TPartIndex::I()->PartIndex(pdg));
  TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(gtrack.fPDG);

  gtrack.SetCharge(part->Charge()/3.);
  gtrack.SetMass(part->Mass());

  if ((HepMC::GenVertex*)genpart->production_vertex())
    {
      gtrack.fXpos = genpart->production_vertex()->position().x();
      gtrack.fYpos = genpart->production_vertex()->position().y();
      gtrack.fZpos = genpart->production_vertex()->position().z();
    }
  else
    {
      gtrack.fXpos = 0;
      gtrack.fYpos = 0;
      gtrack.fZpos = 0;
    }
        
  gtrack.fE = genpart->momentum().e();  //e- 30MeV

  Double_t p = TMath::Sqrt((gtrack.E()-gtrack.Mass())*(gtrack.E()+gtrack.Mass()));
        
  gtrack.SetP(p);
  gtrack.fXdir = genpart->momentum().px()/p;
  gtrack.fYdir = genpart->momentum().py()/p;
  gtrack.fZdir = genpart->momentum().pz()/p;

}


