/*
 * TestConstBzStepper.cxx
 *
 *  Created on: May 8, 2014
 *   Taken from TestConstBzStepper.cxx (written by swenzel)
 */

#include "../inc/ScalarConstFieldHelixStepper.h"
#include "mmalloc.h"
#include "TGeoHelix.h"

#include "GUFieldPropagator.h"

#define ET_TIMER_RESULT_VAR      vx_
#define ET_TIMER_RESULT_VAR_TYPE double*
// #include "Timer.h"

void testGUFieldPropagator(double * px, double * py, double * pz,
                           double * dx, double * dy, double * dz,
                           int * charges,
                           double * mom, double * steps,
                           double * pxn, double * pyn, double * pzn,
                           double * dxn, double * dyn, double * dzn,
                           int np )
{
    GUFieldPropagator stepper(1,1);
    stepper.SetField(0,0,10., false);
    for(int i=0;i<np;++i)
    {
        double const *p;
        double const *d;

        stepper.SetCharge( charges[i] );
        stepper.InitDirection( dx[i], dy[i], dz[i] );
        stepper.InitPoint( px[i], py[i], pz[i] );
        
        const double kB2C_local = -0.299792458e-3;
        // could do a fast square root here
        double dt = sqrt((dx[i]*dx[i]) + (dy[i]*dy[i]));
        double invnorm=1./dt;

        // radius has sign and determines the sense of rotation
        double c = fabs(((kB2C_local*charges[i])*(10.))/mom[i]*dt);
        stepper.SetXYcurvature(c);
        stepper.SetHelixStep(abs(2*M_PI*dz[i])/(c*dt));

        stepper.UpdateHelix();
        stepper.Step( steps[i] );

        d=stepper.GetCurrentDirection();
        p=stepper.GetCurrentPoint();
        pxn[i]=p[0];
        pyn[i]=p[1];
        pzn[i]=p[2];
        dxn[i]=d[0];
        dyn[i]=d[1];
        dzn[i]=d[2];
    }
}

int testPerformance()
{
    int np=128;
    double * px = (double *) _mm_malloc(sizeof(double)*np,32);
    double * py = (double *) _mm_malloc(sizeof(double)*np,32);
    double * pz = (double *) _mm_malloc(sizeof(double)*np,32);
    double * dx = (double *) _mm_malloc(sizeof(double)*np,32);
    double * dy = (double *) _mm_malloc(sizeof(double)*np,32);
    double * dz = (double *) _mm_malloc(sizeof(double)*np,32);

    double * pxn = (double *) _mm_malloc(sizeof(double)*np,32);
    double * pyn = (double *) _mm_malloc(sizeof(double)*np,32);
    double * pzn = (double *) _mm_malloc(sizeof(double)*np,32);
    double * dxn = (double *) _mm_malloc(sizeof(double)*np,32);
    double * dyn = (double *) _mm_malloc(sizeof(double)*np,32);
    double * dzn = (double *) _mm_malloc(sizeof(double)*np,32);

    double * steps = (double *) _mm_malloc(sizeof(double)*np,32);
   // double * B = (double *) _mm_malloc(sizeof(double)*np,32);
    double * mom = (double *) _mm_malloc(sizeof(double)*np,32);
    int * charges = (int *) _mm_malloc(sizeof(int)*np,32);

    // fill data
    for(int i=0;i<np;++i)
    {
      // don't care for units now
      px[i] = 0.;
      py[i] = i;
      pz[i] = i;
      dx[i] = (i%3)? 1 : 0;
      dy[i] = ((i+1)%3)? 1 : 0;
      dz[i] = ((i+2)%3)? 1 : 0;
     // B[i] = 10*i;
      // all doing same step
      steps[i]=10.;
      // momentum
      mom[i]=1234+(i%2);
      charges[i]=(i%2)? 1 : -1;
    }

    double *vx_;
    vx_ = &dxn[0];
    int nreps=1000;

    geantv::ConstBzFieldHelixStepper stepper(10.);
//
//    vx_ = &dyn[0];
//       ET_TIME_THIS( "ROOT performance",
//               testGUFieldPropagator(px,py,pz,dx,dy,dz,charges,mom,steps,
//               pxn,pyn,pzn,dxn,dyn,dzn, np);, nreps, true
//      )
//
//
//    ET_TIME_THIS( "scalar performance",
//            for(int i=0;i<np;++i)
//            {
//               stepper.DoStep(px[i],py[i],pz[i],dx[i],dy[i],dz[i],charges[i],mom[i],steps[i],
//               pxn[i],pyn[i],pzn[i],dxn[i],dyn[i],dzn[i]);
//            }
//    , nreps
//    , true
//    )
//
//    vx_ = &dzn[0];
//    ET_TIME_THIS( "vector performance",
//            stepper.DoStep_v(px,py,pz,dx,dy,dz,charges,mom,steps,
//            pxn,pyn,pzn,dxn,dyn,dzn, np);, nreps, true
//   )

    for(int i=0;i<1000;++i)
    {
       testGUFieldPropagator(px,py,pz,dx,dy,dz,charges,mom,steps,
                             pxn,pyn,pzn,dxn,dyn,dzn, np);
    }

    /***
    for(int i=0;i<1000;++i)
    {
       testTGeoHelix(px,py,pz,dx,dy,dz,charges,mom,steps,
                     pxn,pyn,pzn,dxn,dyn,dzn, np);
    }
    ***/
    
    for(int i=0;i<1000;++i)
    {
       stepper.DoStep_v(px,py,pz,dx,dy,dz,charges,mom,steps,
                        pxn,pyn,pzn,dxn,dyn,dzn, np);
    }
}


int main()
{
    testPerformance();
}

