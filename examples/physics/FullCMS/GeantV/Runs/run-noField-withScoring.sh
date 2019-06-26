###
### GeantV FULL-CMS application script
###
### Note, that the script shows all the possible input arguments that makes possible the configuration
### of the fullCMS GeantV application.
###
###
### ./fullCMS \
$EX_BUILD_DIR/bin/examples/FullCMS/GeantV/FullCMS \
"### detector parameters:"\
  --det-set-gdml                           cms.gdml     "# gdml file "\
"### primary particle generation parameters:"\
  --gun-set-primary-energy                 10           "# kinetic energy of the primary particle [in GeV]"\
  --gun-set-primary-type                   e-           "# primary particle name/type"\
  --gun-set-primary-per-event              10           "# number of primary particles per event"\
"### run configuration parameters:"\
  --config-number-of-buffered-events       10           "# number of events transported at once"\
  --config-total-number-of-events         100           "# total number of events to be simulated"\
  --config-number-of-threads                1           "# number of working threads to be used"\
  --config-number-of-propagators            1           "# number of propagators"\
  --config-run-performance                  0           "# flag to activate performance mode i.e. no scoring"\
  --config-vectorized-geom                  0           "# flag to activate vectorized geometry"\
  --config-external-loop                    0           "# flag to run the application in external loop mode" \
"### basket related parameters:"\
  --config-tracks-per-basket               16           "# default number of tracks per basket" \
  --field-basketized                        0           "# use baskets and vectors for field propoagation " \
"### magnetic field configuration:"\
  --field-type                              0           "# enable magnetic field"

## --gun-set-primary-direction  x=0.1,y=0.9,z=0.1        "# primary particle direction(will be normalized)"\
