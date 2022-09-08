read -p "OK? (y/N): " yes

if [[ ${yes} == "y" ]]; then
    rm -rf log/Submit_preproc_fugaku.* \
       output/* \
       script.*.out \
       script.*.out.*.0 \
       script.*.stats
    rmdir log output
fi
