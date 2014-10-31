//

#ifndef GEANTV_PrimaryGenerator_h
#define GEANTV_PrimaryGenerator_h

class TParticlePDG;
class GeantTrack;

class PrimaryGenerator : public TNamed{
    
public:
    
    // set one GeantTrack primary track properties
    virtual void InitPrimaryGenerator() = 0;
    virtual Int_t NextEvent() = 0;
    virtual void GetTrack(Int_t n, GeantTrack &gtrack) = 0;

};


#endif
