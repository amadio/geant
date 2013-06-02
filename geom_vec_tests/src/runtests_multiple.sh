#!/bin/bash

if [ -n "$1" ]; then
    outdir=$1 
else
    outdir=testresults
fi

compilers="icc g++"
for compiler in $compilers; do
    echo $compiler
    tests=`ls TGeoBBox_Test*SOA*${compiler}`
    echo $tests

    for geotest in $tests; do
	# runtest
	if [ ! -d ${outdir} ]; then
	    mkdir ${outdir}
    	fi

	for run in `seq 0 20`; do
	    outfile="tmp_${geotest}_run${run}.dat"
	    outfile2="out_${geotest}_run${run}.dat"

	    ./${geotest} 2> ${outdir}/${outfile}

	# do some postprocessing / extraction
	    awk '/#/{for (i=2; i<NF; i++){ printf $i"\t"}; print $NF}' ${outdir}/${outfile} > ${outdir}/${outfile2}
	done    
    done
done

