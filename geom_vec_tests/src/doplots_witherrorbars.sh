#!/bin/bash
compilers="icc g++"

if [ -n "$1" ]; then
    outdir=$1 
else
    outdir=testresults_multipleruns
fi

tests=`ls TGeoBBox_Test*SOA*.cxx | xargs -I {} basename {} .cxx`

for geotest in $tests; do
    # produce mean file
    cd $outdir
    for compiler in ${compilers}; do
	python ../statistical_analyse_results.py out_${geotest}_${compiler}_run*.dat 
    done
    cd ..

    gnuplotfile=${geotest}_plot.gpl

    printf "set logscale x\n" > ${outdir}/${gnuplotfile}
    printf "set key spacing 2 top left\n" >> ${outdir}/${gnuplotfile}
    commonopts="w yerrorlines ps 2"
      
    plotcommand="plot"
    for compiler in $compilers; do
	plotfile="out_${geotest}_${compiler}.dat"
	printf "${plotcommand} 'out_${geotest}_${compiler}_mean.dat' u 1:11:12 ${commonopts} pt 5 title '${compiler}_v'\n" >> ${outdir}/${gnuplotfile}
	plotcommand="replot"
	printf "${plotcommand} 'out_${geotest}_${compiler}_mean.dat' u 1:9:10 ${commonopts} pt 4 title '${compiler}_l'\n" >> ${outdir}/${gnuplotfile}

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


