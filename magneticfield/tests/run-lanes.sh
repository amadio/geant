

#  Run multiple runs of testVectorIntegrationDriver with differnt parameters
#
#  Dependencies: 
#     Environment variables: 
#         EX_BUILD_DIR   :  Location of GeantV build
#

runId=$1;  shift ;  ## "runD10"
list=$*

( cd $EX_BUILD_DIR ; make -k -j 4 testVectorIntegrationDriver || ( echo "Make failed - exiting. " ; exit 1 ) )

for lane in $list 
do
   echo "Running with lane = " $lane
   $GVBIN/tests/magneticfield/testVectorIntegrationDriver - ${lane} 2>&1  | tee out-VecIntDrv-ReportLane-no${lane}.${runId}.log | grep ReportOneLane
done
