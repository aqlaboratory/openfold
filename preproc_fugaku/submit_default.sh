NumProcs=1
NumThreads=20

./Submit_preproc_fugaku \
    1:00:00 \
    1 \
    $NumProcs \
    $NumThreads \
    input_examples/4X96.fasta \
    output_$NumThreads \
    uniref90
