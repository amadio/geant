vecLen=$1 ## :512
iter=$2   ## :25000
echo "Running with vector length= $vecLen and $iter iterations"
for i in 1 2 3 4 5 6 7 8 9 10  
do
  ./GUModelBenchmark 256 100000 2>&1 | tee out-vec$vecLen-it$iter.n$i
  sleep 3
done

if [[ ! -f tab-scalar-vec$vecLen-it$iter ]]; then 
   grep Scalar out-vec$vecLen-it$iter.n? > tab-scalar-vec$vecLen-it$iter
   grep Vector out-vec$vecLen-it$iter.n? > tab-vector-vec$vecLen-it$iter
   grep Geant4 out-vec$vecLen-it$iter.n? > tab-alaG4-vec$vecLen-it$iter
fi
