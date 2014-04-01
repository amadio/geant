#ifndef CUDA_GPTRACK_H
#define CUDA_GPTRACK_H 

// interface between G4Track and GPFieldTrack
struct GPTrack { 
  float x; 
  float y;
  float z;
  float px;
  float py;
  float pz;
  float E;
  float q;
  float step;
  float length;
};

#endif
