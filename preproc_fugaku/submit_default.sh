NumThreads=1

./Submit_preproc_fugaku \
    0:10:00 \
    1 \
    $NumThreads \
    input_examples/short_`expr 20 / $NumThreads`.fasta \
    output \
    small_bfd
