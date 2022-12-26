set -eu

# Expected time:
#   small_bfd: 300s
#   mgnify:    1200s
#   uniref90:  1200s
#   pdb70:     300s/#threads

Database=small_bfd
# Database=mgnify
# Database=uniref90
# Database=pdb70

NumNodes=12
NumCpus=48
NumProcs=12
NumThreads=4
./Submit_preproc_fugaku \
    1:00:00 \
    3600 \
    12 \
    $NumProcs \
    $NumThreads \
    /vol0006/share/ra000018/openfold/202212_experiments/pdb_mmcif/latest/pdb_mmcif_sampled_100.fasta \
    output_202212_all \
    $Database \
    ""
