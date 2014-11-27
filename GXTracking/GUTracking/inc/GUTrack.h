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
  int numTracks;
  int *status;
  int *particleType;
  int *id;
  int *parentId;
  int *proc;
  double *x; 
  double *y;
  double *z;
  double *px;
  double *py;
  double *pz;
  double *E;
  double *q;
  double *s;
  double *fPdgV;        // PDG code

} ;

#endif
