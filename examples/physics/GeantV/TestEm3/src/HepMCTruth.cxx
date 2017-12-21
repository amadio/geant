#include "HepMCTruth.h"

#include "GeantTrack.h"

#include "Particle.h"
#include "PhysicsProcess.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/WriterAscii.h"
#ifdef USE_ROOT
#include "HepMC/WriterRoot.h"
#endif
#include "HepMC/Print.h"
#include "Geant/Error.h"

using namespace Geant;


namespace userapplication {

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
  bool HepMCTruth::CheckTrack(GeantTrack &gtrack, MCEvent* /*evt*/)
{
  // check if the track satisfies the condition to be stored
  // 
  bool to_be_stored = false;

  if (gtrack.E() > fEMin) to_be_stored = true; // energy check
  /*
  else if(geantphysics::PhysicsProcess::GetProcessByGlobalIndex(gtrack.Process()))
    {
      if(geantphysics::PhysicsProcess::GetProcessByGlobalIndex(gtrack.Process())->GetName()==
	 "Pair") to_be_stored = true; // process check
    }
  */
  /*
    else if(gtrack.Mother()!=-1 && evt->particles.contains(gtrack.Mother())) // mother energy check
    {
    double motherEnergy = (evt->particles).find(gtrack.Mother())->fE;
      if(motherEnergy > fEMin) to_be_stored = true;
    }
  
  */
  
  return to_be_stored;
}

//______________________________________________________________________________
void HepMCTruth::CloseEvent(int evID) {
  
  HepMC::GenEvent genevt(HepMC::Units::GEV,HepMC::Units::CM);
  genevt.set_event_number(evID);
  
  // map to keep the relation GeantID -> GenParticle
  // needed when connecting daughter to mothers
  std::map<int, HepMC::GenParticlePtr> particle_map;

  HepMC::GenVertexPtr primary_vertex = 0;
  
  auto lt = (events_map.find(evID)->particles).lock_table(); // this lock is probably not needed

  // loop over all particles and create corresponding HepMC::GenParticles
  for (const auto& it : lt) 
    {
      //Printf("particle %i mother %i pdg %i energy %f", it.getKey(), it.getValue()->motherid, it.getValue()->pid, it.getValue()->fE);
      
      HepMC::GenParticlePtr p = HepMC::make_shared<HepMC::GenParticle>(HepMC::FourVector(it.second->fPx,
											 it.second->fPy,
											 it.second->fPz,
											 it.second->fE),
								       geantphysics::Particle::GetParticleByInternalCode(it.second->pid)
								       ->GetPDGCode(),
								       1);
      
      particle_map[it.first] = p;

      // add particle to event
      genevt.add_particle(p);
      
      // check if the particle has end point (vertex) and if yes, create the corresponding HepMC::GenVertex
      // and attach the particle to the vertex
      if(it.second->has_end)
	{	  
	  HepMC::GenVertexPtr vertex = HepMC::make_shared<HepMC::GenVertex>(HepMC::FourVector(it.second->fXend,
											      it.second->fYend,
											      it.second->fZend,
											      it.second->fTend));
	  vertex->add_particle_in(p);
	  genevt.add_vertex(vertex);
	}
    }

  // we need a second loop to connect daughters to mothers
  // we couldn't do it in the previous loop because some of the mothers maybe have been iterated over
  // only after the daughters

  for (const auto& it : lt)
    {
      if(it.second->motherid==-1) // primary particle
	{
	  if(!primary_vertex) // need to create the primary vertex
	    {
	      primary_vertex = HepMC::make_shared<HepMC::GenVertex>(HepMC::FourVector(it.second->fXpos,
										      it.second->fYpos,
										      it.second->fZpos,
										      it.second->fTime));
	      genevt.add_vertex(primary_vertex);
	    }
	  
	  if (it.second->fXpos == primary_vertex->position().x() &&
	      it.second->fYpos == primary_vertex->position().y() &&
	      it.second->fZpos == primary_vertex->position().z())
	    {
	      // coming from 'main' primary vertex
	      primary_vertex->add_particle_out(particle_map[it.first]);      
	    }
	  else
	    {
	      // coming from a 'secondary' primary vertex
	      // need to find/create the corresponding vertex
	      // not handled yet
	      Geant::Printf("HepMCTruth: Looks like you have more than one 'primary' vertex. Handling of such a case is not implemented yet."); 
	    }
	}
      else // not a primary particle
	{
	  // need to find the mother and then attach to its end vertex	  
	  if(particle_map.find(it.second->motherid)!=particle_map.end())
	    {
	      HepMC::GenParticlePtr mother = particle_map[it.second->motherid];
	      
	      if(mother->end_vertex())
		{
		  if( it.second->fXpos == mother->end_vertex()->position().x() &&
		      it.second->fYpos == mother->end_vertex()->position().y() &&
		      it.second->fZpos == mother->end_vertex()->position().z() &&
		      it.second->fTime == mother->end_vertex()->position().t() )		      
		    mother->end_vertex()->add_particle_out(particle_map[it.first]);
		  else
		    {
		      /*
		      Printf("HepMCTruth: Production not at mother's end vertex for PDG %i E %f (skipping particle for the time being)",
			     geantphysics::Particle::GetParticleByInternalCode(it.second->pid)->GetPDGCode(), it.second->fE);
		      */
		      
		      // this is a temporary solution, because the secondary was not produced at the mother end vertex,
		      // but somewhere along the track (like ionization e-)
		      //
		      
		      mother->end_vertex()->add_particle_out(particle_map[it.first]);
		      
		    }
		}
	      else
		{
		  /*
		  Printf("HepMCTruth: Mother %i for particle %i has no end vertex !!!", it.second->motherid,
			 it.first);
		  */
		  
		  // this is a temporary solution
		  // we are attaching the secondaries produced along the track (like ionization e-) to the production vertex
		  // of the mother

		  mother->production_vertex()->add_particle_out(particle_map[it.first]);
		}
	    }
	  else
	    {
	      Printf("HepMCTruth: Not found mother %i for particle %i pdg %i energy %f !!!", it.second->motherid,
		     it.first, geantphysics::Particle::GetParticleByInternalCode(it.second->pid)->GetPDGCode(),
		     it.second->fE); 
	    }
	}
    }
  
  HepMC::Print::listing(genevt);
  //  HepMC::Print::content(genevt);

  output_file->write_event(genevt);
  
  delete events_map.find(evID);
  events_map.erase(evID);  
}

}
