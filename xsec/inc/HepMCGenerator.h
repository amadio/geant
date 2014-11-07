
#ifndef HepMCGenerator_h
#define HepMCGenerator_h

#include "TPartIndex.h"
#include "PrimaryGenerator.h"

#if __cplusplus >= 201103L
#include "HepMC/IO/IO_FileBase.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Search/FindParticles.h"
#endif

class TParticlePDG;
class GeantTrack;


class HepMCGenerator: public PrimaryGenerator{
 private:

#if __cplusplus >= 201103L
    HepMC::IO_Base* input_file;
    HepMC::FindParticles* search;
#endif
    
 public:
    HepMCGenerator();
    HepMCGenerator(std::string& filename);
    ~HepMCGenerator();

  // set one GeantTrack primary track properties
    virtual void InitPrimaryGenerator();
    virtual Int_t NextEvent();
    virtual void GetTrack(Int_t n, GeantTrack &gtrack);

 private:

   HepMCGenerator(const HepMCGenerator &);//no imp.	
   HepMCGenerator& operator=(const HepMCGenerator &);//no imp.

   ClassDef(HepMCGenerator,1)
};

#endif
