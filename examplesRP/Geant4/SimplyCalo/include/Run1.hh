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
/// \file electromagnetic/TestEm3/include/Run11.hh
/// \brief Definition of the Run1 class
//
// $Id: Run1.hh 71375 2013-06-14 07:39:33Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run1_h
#define Run1_h 1

#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "globals.hh"
#include <map>

class DetectorConstruction;
class G4ParticleDefinition;
class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run1 : public G4Run
{
  public:
    Run1(DetectorConstruction*);
   ~Run1();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);

    void FillPerEvent(G4int layerNum, G4int absorNum, G4double edep, G4double stepl);

    virtual void Merge(const G4Run*);
    void         EndOfRun();


    void AddStep();

  private:
    DetectorConstruction*  fDetector;
    G4ParticleDefinition*  fParticle;
    G4double  fEkin;

    G4int fNbOfLayers;
    G4int fNbOfAbsor;

    std::vector<std::vector< G4double> > fSumEdeps;    //[#layerx[#absorbers]]
    std::vector<std::vector< G4double> > fSum2Edeps;   //[#layerx[#absorbers]]
    std::vector<std::vector< G4double> > fSumLengths;  //[#layerx[#absorbers]]
    std::vector<std::vector< G4double> > fSum2Lengths; //[#layerx[#absorbers]]

    unsigned long int   fSumNSteps;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
