#ifndef GUTRACK_H
#define GUTRACK_H 1

struct GUTrack
{
  int status;
  int particleType;
  int id;            // counter
  int parentId;      // id of parent
  int proc;          //  ?? process index ??
  double x;          // x position - rarely relevant
  double y;          // y position - ditto
  double z;          // z position - ditto 
  double px;
  double py;
  double pz;
  double E;
  double q;          // charge ? 
  double s;          // step length ??
  double pdg;        // 
} ;

struct GUTrack_v
{
  int capacity;        // max number of tracks that can be stored
  int numTracks;       // real number of tracks stored
  int *status;         // status of the track: alive or killed (possible at rest ???)
  int *particleType;
  int *id;             
  int *parentId;       // index of the corresponding parent track in GeantTrack_v 
  int *proc;           // process index (not really necessary)
  double *x;           // (x,y,z) position
  double *y;
  double *z;
  double *px;          // momentum (px,py,pz)
  double *py;
  double *pz;
  double *E;           // total energy 
  double *q;           // charge ???
  double *s;           // ???
} ;

#endif
