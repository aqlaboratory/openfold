read -p "OK? (y/N): " yes

if [[ ${yes} == "y" ]]; then
    rm -rf log/Submit_preproc_fugaku.* \
       output/* \
       script.*.out \
       script.*.out.0.* \
       script.*.out.1.* \
       script.*.stats \
       llio_transfer.*.list
    rmdir log output
fi
