# Weak scaling with 1 proc.
NumProcs=1
NumCpus=48
NumThreads=`expr $NumCpus / $NumProcs`

./Submit_preproc_fugaku \
    2:00:00 \
    7200 \
    32 \
    $NumProcs \
    $NumThreads \
    input_examples/short.fasta \
    output_bfd_$NumProcs \
    bfd \
    "--ignore_unique --weak_scale"
