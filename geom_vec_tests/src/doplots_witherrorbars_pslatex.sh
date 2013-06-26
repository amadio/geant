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
#    cd $outdir
#    for compiler in ${compilers}; do
#	python ../statistical_analyse_results.py out_${geotest}_${compiler}_run*.dat 
#    done
#    cd ..

    gnuplotfile=${geotest}_plot.gpl

    printf "set logscale x\n" > ${outdir}/${gnuplotfile}
    printf "set key spacing 3 top left\n" >> ${outdir}/${gnuplotfile}
    printf "set xlabel 'numer of particles'\n" >> ${outdir}/${gnuplotfile}
    printf "set ylabel 'speedup'\n" >> ${outdir}/${gnuplotfile}
    echo "set format '$%g$'" >> ${outdir}/${gnuplotfile}

    commonopts="w yerrorlines ps 2"
      
    plotcommand="plot"
    pt=5
    for compiler in $compilers; do
	plotfile="out_${geotest}_${compiler}.dat"
	printf "${plotcommand} 'out_${geotest}_${compiler}_mean.dat' u 1:11:12 ${commonopts} pt ${pt} title '${compiler}$_$v'\n" >> ${outdir}/${gnuplotfile}
	pt=6
	plotcommand="replot"
#	printf "${plotcommand} 'out_${geotest}_${compiler}_mean.dat' u 1:9:10 ${commonopts} pt 4 title '${compiler}$_$l'\n" >> ${outdir}/${gnuplotfile}
    done
    printf "set term pslatex auxfile color\n"  >> ${outdir}/${gnuplotfile}
    printf "set output 'PerfPlot_${geotest}.tex'\n"  >> ${outdir}/${gnuplotfile}
#    printf "set output 'PerfPlot_${geotest}.ps'\n" 
    printf "replot\n" >> ${outdir}/${gnuplotfile}
#    printf "pause -1" >> ${outdir}/${gnuplotfile}

    # call gnuplot
    cd ${outdir}
    gnuplot ${gnuplotfile}
    # call postprocessing to produce nice latex-like pdf
    fixpslatex_s.pl -i  PerfPlot_${geotest} -o PerfPlotFinal_${geotest} 
    cd $OLDPWD
done


