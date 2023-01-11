set -eu

DataDir=/vol0006/share/ra000018/openfold/202212_experiments/pdb_mmcif/latest
OutDir=output_202212_all

mkdir -p $OutDir

./Submit_preproc \
    128 \
    1:00:00 \
    $DataDir/split \
    $OutDir \
    $DataDir/pdb_mmcif_sampled_10000_mmcif_file_list.txt
