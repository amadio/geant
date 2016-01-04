/*
 * circlestepper.cxx
 *
 *  Created on: May 6, 2014
 *      Author: swenzel
 */

#include <cmath>
#include <iostream>
// #include "vdt/sin.h"
// #include "vdt/cos.h"

inline void steponcircle(int /*charge*/, double R, double x0, double y0, double dx0, double dy0, double step, double &x,
                         double &y, double &dx, double &dy)
{
  double invnorm = 1. / sqrt(dx0 * dx0 + dy0 * dy0);
  double cosa = dx0 * invnorm;
  double sina = dy0 * invnorm;
  double phi = step / R;
  double cosphi = cos(phi);
  double sinphi = sin(phi);

  x = x0 + R * (-sina - (-cosphi * sina - sinphi * cosa));
  y = y0 + R * (cosa - (-sina * sinphi + cosphi * cosa));
  dx = dx0 * cosphi - dy0 * sinphi;
  dy = dx0 * sinphi + dy0 * cosphi;
}

void steponcircle_v(int const *__restrict__ c, double const *__restrict__ R, double const *__restrict__ x0,
                    double const *__restrict__ y0, double const *__restrict__ dx0, double const *__restrict__ dy0,
                    double *__restrict__ step, double *__restrict__ x, double *__restrict__ y, double *__restrict__ dx,
                    double *__restrict__ dy, int np)
{
  for (int i = 0; i < np; ++i) {
    steponcircle(c[i], R[i], x0[i], y0[i], dx0[i], dy0[i], step[i], x[i], y[i], dx[i], dy[i]);
  }
}

void steponhelix(int charge, double R, double x0, double y0, double z0, double dx0, double dy0, double dz0, double step,
                 double &x, double &y, double &z, double &dx, double &dy, double &dz)
{
  // TODO: calculate R and helixstep from physical params
  // helixstep is jump upon a 2pi rotation ( must be a function of B and dz0 )
  double dt = sqrt((dx0 * dx0) + (dy0 * dy0));
  double invnorm = 1. / dt;

  double sina = -charge * dx0 * invnorm;
  double cosa = dy0 * invnorm;
  double helixgradient = dz0 * invnorm * R;

  // can maybe be simplified
  double phi = step / sqrt((R * R) + (helixgradient * helixgradient));

  double xc = x0 - R * cosa;
  double yc = y0 - R * sina;
  double zc = z0;

  double cosphi = cos(phi);
  double sinphi = sin(phi);
  double cosaphi = cosa * cosphi - sina * sinphi;
  double sinaphi = sina * cosphi + cosa * sinphi;

  x = xc + R * (cosaphi);
  y = yc + R * (sinaphi);
  z = zc + helixgradient * phi;
  // let me think: if phi = 2pi
  // then helixstep will come in; so it is
  // something like: z = zc + phi/2pi * Helixstep

  dx = -charge * dt * (sinaphi);
  dy = dt * (cosaphi);
  dz = dz0;
}

int main()
{
  double x0, y0, dx0, dy0, x, y, dx, dy;
  //  double z0=0;
  //  double dz0=0.2;

  x0 = 10;
  y0 = 5;
  dx0 = -0.5;
  dy0 = -0.5;
  double R = 10;

  std::cout << "0"
            << "\t" << x0 << "\t" << y0 << "\t" << dx0 << "\t" << dy0 << "\n";
  for (int i = 1; i <= 60; ++i) {
    steponcircle(1, R, x0, y0, dx0, dy0, 1, x, y, dx, dy);
    std::cout << i << "\t" << x << "\t" << y << "\t" << dx << "\t" << dy << "\n";
    x0 = x;
    y0 = y;
    dx0 = dx;
    dy0 = dy;
  }

  x0 = 10;
  y0 = 5;
  // z0=0.;
  dx0 = -0.5;
  dy0 = -0.5;
  //  dz0 = 0.2;
  std::cout << "0"
            << "\t" << x0 << "\t" << y0 << "\t" << dx0 << "\t" << dy0 << "\n";
  for (int i = 1; i <= 60; ++i) {
    steponcircle(-1, -R, x0, y0, dx0, dy0, 1, x, y, dx, dy);
    std::cout << i << "\t" << x << "\t" << y << "\t" << dx << "\t" << dy << "\n";
    x0 = x;
    y0 = y;
    dx0 = dx;
    dy0 = dy;
  }

  /*
  std::cout << "0" <<"\t" << x0 << "\t" << y0 << "\t" << z0 << " " << dx0 << "\t" << dy0 << "\n";
  for (int i=1;i<=240;++i)
  {
     steponhelix(1, R, x0, y0, z0, dx0, dy0, dz0, 1, x, y, z, dx, dy, dz);
     std::cout << i <<"\t" << x << "\t" << y << "\t" << z << "\t" << dx << "\t" << dy << "\n";
     x0=x;
     y0=y;
     z0=z;
     dx0=dx;
     dy0=dy;
     dz0=dz;
  }
  */
}
