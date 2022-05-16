#!/bin/bash -eux

. scripts/setenv

mkdir -p $WORKDIR

### pytorch
pjsub $PJSUB_OPT_GROUP scripts/install_fugaku_pytorch.sh

### others
pjsub $PJSUB_OPT_GROUP scripts/install_fugaku_others.sh
