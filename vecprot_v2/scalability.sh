#!/bin/bash
# makes scalability plot for up to <maxthr> threads
if [ "$#" -lt 1 ]; then
  echo "Usage: scalability.sh <maxthr> [geomfile] [xsecfile] [fstatefile]"
  exit 1
fi  
if [ "$ROOTSYS" = "" ]; then
 echo "\$ROOTSYS is not defined"
 exit 2
fi 
MAXTHREADS=$1
GEOMFILE="ExN03.root"
XSEC="xsec_FTFP_BERT.root"
FSTATE="fstate_FTFP_BERT.root"
if [ "$2" != "" ]; then
 GEOMFILE=$2
fi
if [ "$3" != "" ]; then
 XSEC=$3
fi
if [ "$4" != "" ]; then
 FSTATE=$4
fi
   
if [ -f scalability.txt ]; then
   rm scalability.txt
fi

echo "using: nthr=$MAXTHREADS geomfile=$GEOMFILE xsec=$XSEC fstate=$FSTATE"
rm -f run*.log

ithr=1
while [ "$ithr" -le "$MAXTHREADS" ]; do
 for count in {1..4} 
 do
  echo -ne "=== threads=$ithr  count=$count ===\r"
  root -b -q run.C\($ithr\,false,\"$GEOMFILE\",\"$XSEC\",\"$FSTATE\"\) >> run$ithr.log 2>>run$ithr.log
 done
 ithr=$(($ithr + 1))
done 

ithr=1
while [ "$ithr" -le "$MAXTHREADS" ]; do
 echo "NTHREADS $ithr" >> scalability.txt
 grep -o 'RT=[.0-9]*' run$ithr.log | sed -r 's/RT=//' >> scalability.txt
 ithr=$(($ithr + 1))
done

$ROOTSYS/bin/root -b -q scalability.C+
exit $?
