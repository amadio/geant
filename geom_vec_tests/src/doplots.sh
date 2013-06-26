#!/bin/bash
compilers="icc g++"

if [ -n "$1" ]; then
    outdir=$1 
else
    outdir=testresults
fi

tests=`ls TGeoBBox_Test*.cxx | xargs -I {} basename {} .cxx`

for geotest in $tests; do
    gnuplotfile=${geotest}_plot.gpl

    printf "set logscale x\n" > ${outdir}/${gnuplotfile}
    commonopts="w lp ps 2"
      
    plotcommand="plot"
    for compiler in $compilers; do
	plotfile="out_${geotest}_${compiler}.dat"
	printf "${plotcommand} 'out_${geotest}_${compiler}.dat' u 1:5 ${commonopts} pt 4 title '${compiler}_l'\n" >> ${outdir}/${gnuplotfile}
	printf "replot 'out_${geotest}_${compiler}.dat' u 1:6 ${commonopts} pt 5 title '${compiler}_v'\n" >> ${outdir}/${gnuplotfile}
	plotcommand="replot"
    done
    printf "set term postscript eps enhanced color\n"  >> ${outdir}/${gnuplotfile}
    printf "set output 'PerfPlot_${geotest}.ps'\n"  >> ${outdir}/${gnuplotfile}
    printf "replot\n" >> ${outdir}/${gnuplotfile}
#    printf "pause -1" >> ${outdir}/${gnuplotfile}

    # call gnuplot
    cd ${outdir}
    gnuplot ${gnuplotfile}
    cd $OLDPWD
done


