#!/bin/bash -eux

. scripts/setenv

mkdir -p $WORKDIR

if [ ! -d $WORKDIR/DeepSpeedCPU ]; then
    echo "Please clone or copy DeepSpeedCPU to $WORKDIR"
fi

### pytorch
pjsub $PJSUB_OPT_GROUP --step --sparam "jnam=install_openfold" scripts/install_fugaku_pytorch.sh

### others
pjsub $PJSUB_OPT_GROUP --step --sparam "jnam=install_openfold,sd=ec!=0:one:0" scripts/install_fugaku_others.sh
