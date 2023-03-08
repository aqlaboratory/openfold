set -eu

#---------- I/O configurations ----------#

# The path to the input .fasta file
InputFile=input_examples/short.fasta

# The path to the output directory
OutputDir=output

# Whether or not LLIO transfer is performed for Python modules, executables, and databases
DoStaging=1

#---------- Job configurations ----------#

# Whether or not each job is submitted
Do_small_bfd=1
Do_mgnify=1
Do_uniref90=1
Do_pdb70=1

# The number of nodes of each search job
NumNodes=1

# The initial number of processes/node and threads/process for each search tool
NumProcsJackhmmer=12
NumThreadsJackhmmer=4
NumProcsHHsearch=4
NumThreadsHHsearch=12

# Limit max memory
LimitMaxMemJackhmmer=1
LimitMaxMemHHsearch=1

# Time limit for each job in the pjsub time format
JobTime_small_bfd=1:00:00
JobTime_mgnify=1:00:00
JobTime_uniref90=1:00:00
JobTime_pdb70=1:00:00

# Time limit for each tool in seconds
Timeout_small_bfd=1800
Timeout_mgnify=3600
Timeout_uniref90=3600
Timeout_pdb70=1800

#----------- Configurations end -----------#


# Submit Jackhmmer jobs
StepNameBase=openfold_preproc_`date +%Y%m%d%H%M%S`
JackhmmerDatabases=(small_bfd mgnify uniref90)
for Database in ${JackhmmerDatabases[@]}; do
    JobTimeVar="JobTime_${Database}"
    TimeoutVar="Timeout_${Database}"
    DoVar="Do_${Database}"
    if (( $DoVar == 1 )); then
	./Submit_preproc_fugaku \
	    ${!JobTimeVar} \
	    ${!TimeoutVar} \
	    $NumNodes \
	    $NumProcsJackhmmer \
	    $NumThreadsJackhmmer \
	    $InputFile \
	    $OutputDir \
	    $Database \
	    "${StepNameBase}_${Database}" \
	    $DoStaging \
	    $LimitMaxMemJackhmmer \
	    ""
    fi
done

# Submit a HHserach job
if (( $Do_pdb70  == 1 )); then
    if (( $Do_uniref90  == 1 )); then
	# Follow the corresponding uniref90 job
	StepName="${StepNameBase}_uniref90"
    else
	# Run independently because the corresponding job is not submitted
	StepName="${StepNameBase}_pdb70"
    fi

    ./Submit_preproc_fugaku \
	$JobTime_pdb70 \
	$Timeout_pdb70 \
	$NumNodes \
	$NumProcsHHsearch \
	$NumThreadsHHsearch \
	$InputFile \
	$OutputDir \
	pdb70 \
	$StepName \
	$DoStaging \
	$LimitMaxMemHHsearch \
	""
fi
