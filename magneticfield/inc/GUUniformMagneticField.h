//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: GUUniformMagField.hh 67605 2013-02-26 20:20:03Z adotti $
//
// 
// class GUUniformMagField
//
// Class description:
//
// Class for creation of Uniform Magnetic Field.

// History:
// - 30.01.97 V.Grichine, Created.
// - 01.08.97 J.Apostolakis, cleanup, new 3-vector constructor, 
//            and removal of helix-stepper (to separate file).
// - 05.11.97 G.Cosmo, added copy constructor and assignment operator.

#ifndef GUUNIFORMMAGFIELD_HH
#define GUUNIFORMMAGFIELD_HH

#include "GUTypes.h"
#include "ThreeVector.h"
#include "GUMagneticField.h"

class GUUniformMagField : public GUMagneticField
{
  public:  // with description
  
    GUUniformMagField(const ThreeVector& FieldVector );
      // A field with value equal to FieldVector.

    GUUniformMagField(double vField,
                      double vTheta,
                      double vPhi     ) ;

    virtual ~GUUniformMagField() ;

    GUUniformMagField(const GUUniformMagField &p);
    GUUniformMagField& operator = (const GUUniformMagField &p);
      // Copy constructor and assignment operator.

    virtual void GetFieldValue(const double yTrack[4],
                                     double *MagField) const ;

    void SetFieldValue(const ThreeVector& newFieldValue);

    ThreeVector GetConstantFieldValue() const;
      // Return the field value
    
    virtual GUUniformMagField* Clone() const;

  private:

    double fFieldComponents[3] ;
};

#endif
