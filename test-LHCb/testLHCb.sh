#!/bin/bash

export VP_DATA_XSEC=$GEANTV_DATA_DIR/xsec_FTFP_BERT_G496p02_1mev.root 
export VP_DATA_FSTA=$GEANTV_DATA_DIR/fstate_FTFP_BERT_G496p02_1mev.root

#################################################################
# GDML Geomtery File 
GEOMFILE="$GEANTV_DATA_DIR/LHCb_201603.gdml"
#################################################################

#################################################################
# ROOT file that contains the pre-generated primary events
EVENTFILE="$GEANTV_DATA_DIR/pp14TeVminbias.root"
#################################################################

#################################################################
# Geant4 macro file with some Geant4 commands
GEANTMACRO="g4macro_TPHYS.mac" 
#################################################################

#################################################################
# Low energy cut value given in GeV units
LOWENERGYCUT="0.001"
#################################################################

#################################################################
# (OPTIONAL) Physics list name: Tabulated physics is the default
PHYSLISTNAME="TABPHYS"
#################################################################

#################################################################
# Set the required score type (0 is default i.e. no scoring)
SCORETYPEFLAG="2"
#################################################################
 
../../installation/bin/testlhcb \
--geomFile $GEOMFILE \
--eventFile $EVENTFILE \
--geantMacro $GEANTMACRO \
--lowEnergyCut $LOWENERGYCUT \
--physListName $PHYSLISTNAME \
--scoreType $SCORETYPEFLAG

