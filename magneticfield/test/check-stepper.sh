
#  Simple exercise script to check results for driver (chosen Cash Karp.)
#  
#  Created:   J. Apostolakis,  20 October 2017
#

dat=`date "+%Y-%m-%d-%H:%M:%S"`

outf="out-stepperCK-B3-s250-n17" 
outNew=${outf}.new

if [[ -f $outNew ]]; then 
   echo "Output file $outNew exists - deleting it. " 
   sleep 2
   rm  -f $outNew
fi

##                             Driver Step  #steps B-field
$GVBIN/testStepperFixedCashKarp 5  250.0   17     3  > $outNew

diff $outf  $outNew   # ${outf}.new
mv ${outf}.new ${outf}.$dat

##  Run with larger step size, a multiple of the previous one (eventual comparison)
##
outf2="out-stepperCK-B3-s2k-n6" 
outNew2="${outf2}.new"
if [[ -f $outNew2 ]]; then 
   echo "Output file $outNew2 exists - deleting it. " 
   sleep 2
   rm  -f $outNew2
fi

##                             Driver Step  #steps B-field
$GVBIN/testStepperFixedCashKarp 5 2000.0    6     3  > $outNew2   ## out-stepperCK-B3-s2k-n6

diff $outf2 $outNew2   # ${outf2}.new
mv ${outf2}.new ${outf2}.$dat
