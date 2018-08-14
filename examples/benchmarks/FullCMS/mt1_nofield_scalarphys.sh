#!/bin/bash
### gperftools file name
export GEANT_PERFTOOLS_FILE="prof_mt1_nofield_scalarphys.out"

### Primary generator
GUN_PRIMARY_TYPE="e-"
GUN_PRIMARY_ENERGY=100 #GeV

### Field
FIELD_ACTIVE=0
FIELD_USE_RK=0

### Configuration
CONFIG_THREADS=1
CONFIG_EVENTS=10
CONFIG_PRIMARIES=10
CONFIG_PERFORMANCE=1
CONFIG_BASKETIZED_FIELD=0
CONFIG_BASKETIZED_PHYSICS=0
CONFIG_BASKETIZED_GEOM=0
CONFIG_BASKETIZED_MSC=0
