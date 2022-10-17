/home/acd13359ld/venv/openfold_fugaku/tools_install/bin/hhblits \
    -i \
    /home/acd13359ld/repo/openfold_fugaku/preproc_fugaku/input_examples/short.fasta \
    -cpu \
    20 \
    -oa3m \
    output.a3m \
    -n \
    3 \
    -e \
    0.001 \
    -maxseq \
    1000000 \
    -realign_max \
    100000 \
    -maxfilt \
    100000 \
    -min_prefilter_hits \
    1000 \
    -d \
    /groups/4/gcd50734/openfold/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    -d \
    /groups/4/gcd50734/openfold/uniclust30/uniclust30_2018_08/uniclust30_2018_08


# -o \
    # /dev/null \
