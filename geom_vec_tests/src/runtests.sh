#!/bin/bash

if [ -n "$1" ]; then
    outdir=$1 
else
    outdir=testresults
fi

compilers="icc g++"
for compiler in $compilers; do
    echo $compiler
    tests=`ls TGeoBBox_Test*${compiler}`
    echo $tests

    for geotest in $tests; do
	# runtest
	outfile="tmp_${geotest}.dat"
	outfile2="out_${geotest}.dat"
	if [ ! -d ${outdir} ]; then
	    mkdir ${outdir}
    	fi
	./${geotest} 2> ${outdir}/${outfile}

	# do some postprocessing / extraction
	awk '/#/{for (i=2; i<NF; i++){ printf $i"\t"}; print $NF}' ${outdir}/${outfile} > ${outdir}/${outfile2}
    done
done

