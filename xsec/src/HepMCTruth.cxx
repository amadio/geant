#include "HepMCTruth.h"

#include "GeantTrack.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/WriterAscii.h"
#ifdef USE_ROOT
#include "HepMC/WriterRoot.h"
#endif
#include "HepMC/Print.h"
#include "Geant/Error.h"

using namespace Geant;

#ifdef USE_ROOT
ClassImp(HepMCTruth)
#endif

//______________________________________________________________________________
HepMCTruth::HepMCTruth() : output_file(0) 
{
}

//______________________________________________________________________________
HepMCTruth::HepMCTruth(std::string &filename) : output_file(0), fEMin(0) {
  if (filename.substr(filename.find_last_of(".") + 1) == "hepmc3") {
    output_file = new HepMC::WriterAscii(filename);
  }
#ifdef USE_ROOT
  else if (filename.substr(filename.find_last_of(".") + 1) == "root") {
    output_file = new HepMC::WriterRoot(filename);
  }
#endif
  else {
    std::cout << "Unrecognized filename extension (must be .hepmc3 or .root)" << std::endl;
  }
}

//______________________________________________________________________________
HepMCTruth::~HepMCTruth() {
  output_file->close();
  delete output_file;
}

//______________________________________________________________________________
void HepMCTruth::InitMCTruthMgr() {}

//______________________________________________________________________________
bool HepMCTruth::CheckTrack(GeantTrack &gtrack, MCEvent* evt)
{
  // check if the track satisfies the condition to be stored
  // 
  bool to_be_stored = false;
  
  if (gtrack.E() > fEMin) to_be_stored = true;
  else if(gtrack.Mother()!=0 && evt->particles.contains(gtrack.Mother()))
    {
      double motherEnergy = (evt->particles).find(gtrack.Mother())->fE;
      if(motherEnergy > fEMin) to_be_stored = true;
    }
  return to_be_stored;
}

//______________________________________________________________________________
void HepMCTruth::CloseEvent(int evID) {
  
  HepMC::GenEvent genevt(HepMC::Units::GEV,HepMC::Units::MM);
  genevt.set_event_number(evID);
  
  // map to keep the relation GeantID -> GenParticle
  // needed when connecting daughter to mothers
  std::map<int, HepMC::GenParticlePtr> particle_map;
  
  HepMC::GenVertexPtr primary_vertex = HepMC::make_shared<HepMC::GenVertex>();
  genevt.add_vertex(primary_vertex);
 
  auto lt = (events_map.find(evID)->particles).lock_table();
  for (const auto& it : lt) 
    {
      //      Printf("particle %i mother %i pdg %i energy %f", it.getKey(), it.getValue()->motherid, it.getValue()->pid, it.getValue()->fE);
      
      HepMC::GenParticlePtr p = HepMC::make_shared<HepMC::GenParticle>(HepMC::FourVector(it.second->fPx,
											 it.second->fPy,
											 it.second->fPz,
											 it.second->fE),
								       it.second->pid, 1);

      particle_map[it.first] = p;
		   
      // check if the particle has end point (vertex)
      if(it.second->has_end)
	{	  
	  HepMC::GenVertexPtr vertex = HepMC::make_shared<HepMC::GenVertex>(HepMC::FourVector(it.second->fXend,
											      it.second->fYend,
											      it.second->fZend,
											      it.second->fTend));
	  vertex->add_particle_in(p);
	  genevt.add_vertex(vertex);
	}
      
      // add particle to event
      genevt.add_particle(p);
      
      // now look at the production vertex...
      if(it.second->motherid==0)
	{
	  // primary particle
	  if (it.second->fXpos==0 && it.second->fYpos==0 && it.second->fZpos==0)
	    {
	      // coming from 'main' primary vertex
	      primary_vertex->add_particle_out(p);      
	    }
	  else
	    {
	      // coming form 'secondary' primary vertex
	      // need to find/create the corresponding vertex

	      Geant::Printf("HepMCTruth: Primary particle not coming from (0,0,0) <- to be done");

	      // check if the vertex already exists
	    }
	}
      else
	{
	  // secondary particle
	  // need to find the mother and then somehow attach to its end vertex
	  
	  if(particle_map.find(it.second->motherid)!=particle_map.end())
	    {
	      HepMC::GenParticlePtr mother = particle_map[it.second->motherid];
	  
	      if(mother->end_vertex())
		{
		  if( it.second->fXpos == mother->end_vertex()->position().x() &&
		      it.second->fYpos == mother->end_vertex()->position().y() &&
		      it.second->fZpos == mother->end_vertex()->position().z() &&
		      it.second->fTime == mother->end_vertex()->position().t() )		      
		    mother->end_vertex()->add_particle_out(p);
		  else
		    {
		      //Printf("HepMCTruth: Production not at mother's end vertex for PDG %i E %f (skipping particle for the time being)",
		      //     it.second->pid, it.second->fE);
		    }
		}
	      else {
		//Printf("HepMCTruth: Mother %i for particle %i has no end vertex !!!", it.second->motherid,
		//	  it.first);
	      }
	    }
	  else {
	    //Printf("HepMCTruth: Not found mother %i for particle %i pdg %i energy %f !!!", it.second->motherid,
	    //      it.first,
	    //	      it.second->pid, it.second->fE); 
	  }
	}
    }  
  HepMC::Print::listing(genevt);
  //  HepMC::Print::content(genevt);

  output_file->write_event(genevt);
  
  delete events_map.find(evID);
  events_map.erase(evID);  
}
