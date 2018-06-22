GeantPhysicsTest(
  INCDIRS physics/neutron/nudy/inc physics/hadron/cross_sections/inc physics/hadron/models/inc physics/hadron/processes/inc
          physics/neutron/models/inc physics/neutron/xsec/inc
)
# or could do
#   GeantPhysicsTest(MAIN test_stem SCRIPTS README)
# if we have test_stem.cc and this will produce a exec named "test_stem"
#  INCDIRS physics/hadron/cross_sections/inc physics/hadron/models/inc physics/hadron/processes/inc
