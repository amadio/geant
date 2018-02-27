#ifndef MSC_DATA__H
#define MSC_DATA__H

namespace geantphysics {

struct MSCdata {
  double fLambda0 = 0;                  /* elastic mean free path */
  double fLambda1 = 0.;                 /* first transport mean free path */
  double fScrA = 0.;                    /* screening parameter if any */
  double fG1 = 0.;                      /* first transport coef. */
  double fRange = 1.e+21;               /* range of the particle */

  double fTheInitialRange = 0.;         /* the initial (first step or first step in volume) range value of the particle */
//  double fTheRangeFactor = 0.;          /* a step limit factor set */
  double fTheTrueStepLenght = 0.;       /* the true step length */
  double fTheTransportDistance = 0.;    /* the straight line distance between the pre- and true post-step points */
  double fTheZPathLenght = 0.;          /* projection of transport distance along the original direction */
  double fTheTrueGeomLimit = 1.e+20;    /* geometrical limit converted to true step length */
  double fTheDisplacementVectorX = 0.;  /* displacement vector X component */
  double fTheDisplacementVectorY = 0.;  /* displacement vector Y component */
  double fTheDisplacementVectorZ = 0.;  /* displacement vector Z component */
  double fTheNewDirectionX = 0.;        /* new direction X component */
  double fTheNewDirectionY = 0.;        /* new direction Y component */
  double fTheNewDirectionZ = 1.;        /* new direction Z component */
  double fPar1 = -1.;                   /* ??? */
  double fPar2 = 0.;                    /* ??? */
  double fPar3 = 0.;                    /* ??? */

  // Status flags
  bool fIsEverythingWasDone = false;    /* indicates if everything could be done in the step limit phase */
  bool fIsMultipleSacettring = false;   /* indicates that msc needs to be perform (i.e. compute angular deflection) */
  bool fIsSingleScattering = false;     /* indicates that single scattering needs to be done */
  bool fIsEndedUpOnBoundary = false;    /* indicates that geometry was the winer */
  bool fIsNoScatteringInMSC = false;    /* indicates that no scattering happend (i.e. forward) in msc */
  bool fIsNoDisplace = false;           /* indicates that displacement is not computed */
  bool fIsInsideSkin = false;           /* indicates that the particle is within skin from/to boundary */
  bool fIsWasOnBoundary = false;        /* indicates that boundary crossing happend recently */
  bool fIsFirstStep = true;             /* indicates that the first step is made with the particle */
  bool fIsFirstRealStep = false;        /* indicates that the particle is making the first real step in the volume i.e. */

  /**
   * @brief Function that set X, Y, Z component of the displacement vector provided by msc
   *
   * @param x X position
   * @param y Y position
   * @param z Z position
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetDisplacement(double x, double y, double z) {
    fTheDisplacementVectorX = x;
    fTheDisplacementVectorY = y;
    fTheDisplacementVectorZ = z;
  }

  /**
   * @brief Function that set the new X, Y, Z directions proposed by msc (will be applied or not depending other conditions)
   *
   * @param dx X direction
   * @param dy Y direction
   * @param dz Z direction
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNewDirectionMsc(double dx, double dy, double dz) {
    fTheNewDirectionX = dx;
    fTheNewDirectionY = dy;
    fTheNewDirectionZ = dz;
  }

  void Print() const {
    printf("fLambda0=%g fLambda1=%g fScrA=%g fG1=%g fRange=%g fTheInitialRange=%g \
fTheTrueStepLenght=%g fTheTransportDistance=%g fTheZPathLenght=%g fTheTrueGeomLimit=%g \
fTheDisplacementVectorX=%g fTheDisplacementVectorY=%g fTheDisplacementVectorZ=%g \
fTheNewDirectionX=%g fTheNewDirectionY=%g fTheNewDirectionZ=%g fPar1=%g fPar2=%g fPar3=%g \
fIsEverythingWasDone=%d fIsMultipleSacettring=%d fIsSingleScattering=%d fIsEndedUpOnBoundary=%d \
fIsNoScatteringInMSC=%d fIsNoDisplace=%d fIsInsideSkin=%d fIsWasOnBoundary=%d \
fIsFirstStep=%d fIsFirstRealStep=%d\n", fLambda0, fLambda1, fScrA, fG1, fRange, fTheInitialRange,
           fTheTrueStepLenght, fTheTransportDistance, fTheZPathLenght, fTheTrueGeomLimit,
           fTheDisplacementVectorX, fTheDisplacementVectorY, fTheDisplacementVectorZ,
           fTheNewDirectionX, fTheNewDirectionY, fTheNewDirectionZ, fPar1, fPar2, fPar3,
           fIsEverythingWasDone, fIsMultipleSacettring, fIsSingleScattering, fIsEndedUpOnBoundary,
           fIsNoScatteringInMSC, fIsNoDisplace, fIsInsideSkin, fIsWasOnBoundary,
           fIsFirstStep, fIsFirstRealStep);
  }
};

}

#endif  // MSC_DATA__H
