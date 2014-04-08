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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronBuilder_WP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 13.06.2006 G.Folger: (re)move elastic scatterring 
//
//----------------------------------------------------------------------------
//
#include "G4NeutronBuilder_WP.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "TabulatedHadronProcess.hh"

G4NeutronBuilder_WP::
G4NeutronBuilder_WP(): wasActivated(false) 
{
  theNeutronInelastic = new G4NeutronInelasticProcess;
  theNeutronCapture = new G4HadronCaptureProcess;
  theNeutronFission = new G4HadronFissionProcess;

  theWrappedNeutronInelastic = 
    new TabulatedHadronProcess(theNeutronInelastic->GetProcessName(),
                               theNeutronInelastic->GetProcessType());
  theWrappedNeutronCapture   = 
    new TabulatedHadronProcess(theNeutronCapture->GetProcessName(),
                               theNeutronCapture->GetProcessType());
  theWrappedNeutronFission   = 
    new TabulatedHadronProcess(theNeutronFission->GetProcessName(),
                               theNeutronFission->GetProcessType());
}

G4NeutronBuilder_WP::
~G4NeutronBuilder_WP() 
{
  delete theNeutronInelastic;
  delete theNeutronCapture;
  delete theNeutronFission;

  delete theWrappedNeutronInelastic;
  delete theWrappedNeutronCapture;
  delete theWrappedNeutronFission;
}

void G4NeutronBuilder_WP::
Build()
{
  wasActivated = true;
  std::vector<G4VNeutronBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(theNeutronInelastic);
    (*i)->Build(theNeutronCapture);
    (*i)->Build(theNeutronFission);
  }
  G4ProcessManager * theProcMan = G4Neutron::Neutron()->GetProcessManager();
  //  theProcMan->AddDiscreteProcess(theNeutronInelastic);
  theWrappedNeutronInelastic->RegisterProcess(theNeutronInelastic);
  theProcMan->AddDiscreteProcess(theWrappedNeutronInelastic);

  //  theProcMan->AddDiscreteProcess(theNeutronCapture);
  theWrappedNeutronCapture->RegisterProcess(theNeutronCapture);
  theProcMan->AddDiscreteProcess(theWrappedNeutronCapture);

  //  theProcMan->AddDiscreteProcess(theNeutronFission);
  theWrappedNeutronFission->RegisterProcess(theNeutronFission);
  theProcMan->AddDiscreteProcess(theWrappedNeutronFission);
}
// 2002 by J.P. Wellisch
