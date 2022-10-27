ToolTimeLimit=3600

if "true"; then
    # FIXME: tool time limit
    ./Submit_preproc_fugaku \
	1:00:00 \
	$ToolTimeLimit \
	1 \
	1 \
	20 \
	input_examples/short.fasta \
	output \
	small_bfd \
	""

elif "true"; then
    dbs=(uniref90 mgnify small_bfd bfd)
    for db in ${dbs[@]}; do
	./Submit_preproc_fugaku \
	    1:00:00 \
	    $ToolTimeLimit \
	    1 \
	    1 \
	    20 \
	    input_examples/short.fasta \
	    output \
	    ${db} \
	    ""
	sleep 1
    done

elif "true"; then
    # Strong scaling with 1 proc.
    NumProcs=1
    NumThreadss=(1 2 5 10 20) # ABCI
    # NumThreadss=(1 2 6 12 24 48) # Fugaku
    for NumThreads in ${NumThreadss[@]}; do
	./Submit_preproc_fugaku \
	    1:00:00 \
	    $ToolTimeLimit \
	    1 \
	    $NumProcs \
	    $NumThreads \
	    input_examples/short.fasta \
	    output_strong_$NumThreads \
	    small_bfd \
	    "--ignore_unique"
	sleep 1
    done

else
    # Weak scaling with 1 proc.
    source ../scripts/setenv # required for getting NUM_CORES
    NumProcss=(1 2 5 10 20) # ABCI
    # NumProcss=(1 2 6 12 24 48) # Fugaku
    for NumProcs in ${NumProcss[@]}; do
	NumCpus=$NUM_CORES
	NumThreads=`expr $NumCpus / $NumProcs`

	if [[ ! `expr $NumProcs \* $NumThreads` == $NumCpus ]]; then
	    echo "NumCpus ($NumCpus) must be dividible by NumProcs ($NumProcs) and NumThreads ($NumThreads)" >&2
	    exit
	fi

	./Submit_preproc_fugaku \
	    1:00:00 \
	    $ToolTimeLimit \
	    1 \
	    $NumProcs \
	    $NumThreads \
	    input_examples/short.fasta \
	    output_weak_$NumProcs \
	    small_bfd \
	    "--ignore_unique --weak_scale"
	sleep 1
    done
fi
