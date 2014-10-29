#ifndef GUTRACK_H
#define GUTRACK_H 1

struct GUTrack
{
  int status;
  int id;
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
  int *status;
  int *id;
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
