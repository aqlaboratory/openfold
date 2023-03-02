set -eu

StepNameBase=openfold_preproc_`date +%Y%m%d%H%M%S`
Databases=(small_bfd mgnify uniref90)
DoStaging=0

for Database in ${Databases[@]}; do
    ./Submit_preproc_fugaku \
	1:00:00 \
	3600 \
	1 \
	12 \
	4 \
	input_examples/short.fasta \
	output \
	$Database \
	${StepNameBase}_${Database} \
	$DoStaging \
	""
done

./Submit_preproc_fugaku \
    1:00:00 \
    3600 \
    1 \
    12 \
    4 \
    input_examples/short.fasta \
    output \
    pdb70 \
    ${StepNameBase}_uniref90 \
    $DoStaging \
    ""
