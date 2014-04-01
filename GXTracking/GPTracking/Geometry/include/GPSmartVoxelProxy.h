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
// $Id: G4SmartVoxelProxy.hh,v 1.8 2006-06-29 18:32:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4SmartVoxelProxy
//
// Class description:
//
// Class for proxying smart voxels. The class represents either a header
// (in turn refering to more VoxelProxies) or a node. If created as a node,
// calls to GetHeader cause an exception, and likewise GetNode when a header.
//
// Note that the proxy does NOT gain deletion responsibility for proxied
// objects.

// History:
// 12.07.95 P.Kent Initial version
// 03.08.95 P.Kent Updated to become non abstract class, removing
//                 HeaderProxy and NodeProxy derived classes
// --------------------------------------------------------------------
#ifndef GPSMARTVOXELPROXY_HH
#define GPSMARTVOXELPROXY_HH

#include "GPTypeDef.h"
//#include <assert.h>

#include "GPSmartVoxelNode.h"
#include "GPSmartVoxelHeader.h"

struct GPSmartVoxelProxy 
{
  GEOMETRYLOC struct GPSmartVoxelHeader* fHeader;
  GEOMETRYLOC struct GPSmartVoxelNode* fNode;
};

extern "C" {

FQUALIFIER
G4bool GPSmartVoxelProxy_IsHeader(GEOMETRYLOC GPSmartVoxelProxy *This);

FQUALIFIER
G4bool GPSmartVoxelProxy_IsNode(GEOMETRYLOC GPSmartVoxelProxy *This);

FQUALIFIER
GEOMETRYLOC GPSmartVoxelNode* GPSmartVoxelProxy_GetNode(GEOMETRYLOC GPSmartVoxelProxy *This);

FQUALIFIER
GEOMETRYLOC GPSmartVoxelHeader* GPSmartVoxelProxy_GetHeader(GEOMETRYLOC GPSmartVoxelProxy *This);

FQUALIFIER
G4bool GPSmartVoxelProxy_operator_equal(GPSmartVoxelProxy *This,
                                        GPSmartVoxelProxy *v); 

}
#endif
