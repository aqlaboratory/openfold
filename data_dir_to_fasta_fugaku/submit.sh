pjsub -L node=1 \
      -L elapse=0:30:00 \
      -L rscgrp=small \
      -L freq=2200,eco_state=2 \
      ${PJSUB_OPT} \
      -j \
      ./job.sh
