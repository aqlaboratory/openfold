# Expected time:
#   small_bfd: 300s
#   mgnify:    1200s
#   uniref90:  1200s
#   pdb70:     300s/#threads

Database=small_bfd
# Database=mgnify
# Database=uniref90
# Database=pdb70

NumCpus=48
NumProcss=(1 2 6 12 24 48)
for NumProcs in ${NumProcss[@]}; do
    NumThreads=`expr $NumCpus / $NumProcs`
    ./Submit_preproc_fugaku \
	1:00:00 \
	3600 \
	1 \
	$NumProcs \
	$NumThreads \
	input_examples/202212_experiments/20221219_pdb_mmcif_sampled_10000_p50.fasta \
	output_202212_thread_tuning_p$NumProcs \
	$Database \
	"--ignore_unique --weak_scale"
    sleep 1
done
