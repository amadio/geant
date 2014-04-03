#ifndef SOA_KERNEL_H
#define SOA_KERNEL_H

//typedef float G4double;
typedef double G4double;

struct StructSoA4 {
  G4double *x;
  G4double *y;
  G4double *z;
  G4double *t;
};

struct StructAoS4 {
  G4double  x;
  G4double  y;
  G4double  z;
  G4double  t;
};

struct StructSoA8 {
  G4double *x;
  G4double *y;
  G4double *z;
  G4double *t;
  G4double *u;
  G4double *v;
  G4double *w;
  G4double *s;
};

struct StructAoS8 {
  G4double  x;
  G4double  y;
  G4double  z;
  G4double  t;
  G4double  u;
  G4double  v;
  G4double  w;
  G4double  s;
};

void soa4_gpu(StructSoA4  soa, size_t numTracks, int NB, int NT);
void soa4_cpu(StructSoA4  soa, size_t numTracks);

void aos4_gpu(StructAoS4 *aos, size_t numTracks, int NB, int NT);
void aos4_cpu(StructAoS4 *aos, size_t numTracks);

void soa8_gpu(StructSoA8  soa, size_t numTracks, int NB, int NT);
void soa8_cpu(StructSoA8  soa, size_t numTracks);

void aos8_gpu(StructAoS8 *aos, size_t numTracks, int NB, int NT);
void aos8_cpu(StructAoS8 *aos, size_t numTracks);

#endif
