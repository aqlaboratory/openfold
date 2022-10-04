NumProcs=1

# Strong scaling with 1 proc.
NumThreadss=(1 2 5 10 20)
for NumThreads in ${NumThreadss[@]}; do
    ./Submit_preproc_fugaku \
	1:00:00 \
	1 \
	$NumProcs \
	$NumThreads \
	input_examples/short/short.fasta \
	output_$NumThreads \
	uniref90
    sleep 1
done
