read -p "OK? (y/N): " yes

if [[ ${yes} == "y" ]]; then
    rm -rf \
       job.sh.*.out \
       llio_transfer.*.list
fi
