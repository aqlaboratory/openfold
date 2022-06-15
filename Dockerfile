FROM nvidia/cuda:10.2-cudnn8-runtime-ubuntu18.04

RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub

RUN apt-get update \
    && apt-get install -y wget cuda-minimal-build-10-2 git zip \
    && apt-get clean

RUN wget -q -P /tmp \
  https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
  && bash /tmp/Mambaforge-Linux-x86_64.sh -b -p /opt/mamba \
  && rm /tmp/Mambaforge-Linux-x86_64.sh

ENV PATH="/opt/mamba/bin:$PATH"

RUN wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm awscliv2.zip

RUN git clone --branch main --single-branch https://github.com/aqlaboratory/openfold.git /opt/openfold \
    && cd /opt/openfold \
    && git reset ec5619fc970e28e7b81ce452f5e08e7dd6a7cb31 \
    && rm -rf /opt/openfold/imgs /opt/openfold/notebooks /opt/openfold/tests

RUN mamba env update -n base --file /opt/openfold/environment.yml \
    && mamba clean --all

RUN wget -q -P /opt/openfold/openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

RUN patch -p0 -d /opt/mamba/lib/python3.7/site-packages/ < /opt/openfold/lib/openmm.patch

COPY run_batch_job.sh /opt/openfold

WORKDIR /opt/openfold 

RUN python3 setup.py install