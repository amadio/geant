#ifndef GUTRACK_H
#define GUTRACK_H 1

struct GUTrack
{
  int status;
  int particleType;
  int id;            // counter
  int parentId;      // id of parent
  int proc;
  double x; 
  double y;
  double z;
  double px;
  double py;
  double pz;
  double E;
  double q;
  double s;
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
} ;

#endif
