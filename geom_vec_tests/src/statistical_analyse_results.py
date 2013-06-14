"""
A little script that calculates means and errors from benchmark results
coming from different runs
Basically this is just implementing a mean / error operation on matrices
using numpy
"""
import numpy
import sys
import math
import re

matrixlist=list()
filenamelist=list()

def fill_matrices():
    """ 
    reads in result data from a file and stores it in a numpy matrix
    the reference to the matrix will be kept in the matrixlist
    """
    for datafile in filenamelist:
        print datafile
        M=numpy.loadtxt(datafile)
        matrixlist.append(M)


def determine_dimensions():
    """
    a simple consistency check on all the data
    """
    None


def getmean_and_error(i, j):
    """
    calculates mean and error of the data in row i, column j
    returning a tuple
    """
    mean=0.;
    meansqr=0.;
    numberm=len(matrixlist)
    for m in matrixlist:
        mean = mean + m[i][j]
        meansqr = meansqr + m[i][j]*m[i][j]
    mean=mean/numberm;
    meansqr=meansqr/numberm;
    return (mean, math.sqrt((meansqr - mean*mean)))


def getoutfilename():
    """ 
    generally we assume that all filename have same fileroot and a suffix run%d.dat
    TODO: give command line option with -o ...
    returns: fileroot_mean.dat
    """
    # get filename from first file in filelist
    filename=filenamelist[0]
    matcher=re.compile('(.*_)(run[0-9]+)\.dat$')
    result=matcher.match(filename)
    if result:
        root=result.group(1);
        return root+"mean.dat"
    return ""


def iterate():
    rownumber=len(matrixlist[0])
    columnumber=len(matrixlist[0][0])
   
    output=numpy.ndarray( (rownumber, 2*columnumber ) )
    for i in range(rownumber):
        for j in range(columnumber):
            (mean,error)=getmean_and_error(i,j)
            output[i][2*j]=mean
            output[i][2*j+1]=error            
                       
    return output

def main():
    
    # we are expecting all relevant result files as input
    for f in sys.argv[1:]:
        filenamelist.append(f)

    fill_matrices()

    print matrixlist
    print getmean_and_error(3,3)

    outputmatrix=iterate()
    numpy.savetxt(getoutfilename(),outputmatrix,'%.6g')                      

if __name__ == "__main__":
    main()

