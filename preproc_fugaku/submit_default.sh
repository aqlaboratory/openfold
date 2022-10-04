# # Strong scaling with 1 proc.
# NumProcs=1
# NumThreadss=(1 2 5 10 20)
# for NumThreads in ${NumThreadss[@]}; do
#     ./Submit_preproc_fugaku \
    #	1:00:00 \
    #	1 \
    #	$NumProcs \
    #	$NumThreads \
    #	input_examples/short/short.fasta \
    #	output_$NumThreads \
    #	uniref90
#     sleep 1
# done

# Weak scaling with 1 proc.
NumProcss=(1 2 5 10 20)
for NumProcs in ${NumProcss[@]}; do
    NumCpus=20
    NumThreads=`expr $NumCpus / $NumProcs`

    if [[ ! `expr $NumProcs \* $NumThreads` == $NumCpus ]]; then
	echo "NumCpus ($NumCpus) must be dividible by NumProcs ($NumProcs) and NumThreads ($NumThreads)" >&2
	exit
    fi

    ./Submit_preproc_fugaku \
	1:00:00 \
	1 \
	$NumProcs \
	$NumThreads \
	input_examples/4X96.fasta \
	output_$NumProcs \
	uniref90
    sleep 1
done
