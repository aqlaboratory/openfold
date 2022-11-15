NumProcs=1
NumCpus=40
NumThreads=`expr $NumCpus / $NumProcs`

./Submit_preproc_fugaku \
    1:00:00 \
    3600 \
    1 \
    $NumProcs \
    $NumThreads \
    input_examples/short.fasta \
    output_bfd_$NumProcs \
    bfd_deflate \
    ""
