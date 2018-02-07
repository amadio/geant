   ln -s ../../VecMagFieldRoutine/cmsmagfield2015.txt .
   $GVBIN/testMagFieldEquation
   $GVBIN/testConstVecFieldStepper  1 10 1.0 0.0 0.0  | egrep -v '(Bad|bad)'

   $GVBIN/testStepperFixed
   $GVBIN/testStepperFixedCashKarp
   $GVBIN/testIntegrationDriver > newout
   diff out-testIntegrationDriver.def newout

$GVBIN/testVectorMagFieldEquation || echo "Failed test: testVectorMagFieldEquation  ( Equation for 'vector' of tracks in Magnetic Field ) "
