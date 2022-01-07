FROM nvidia/cuda:11.3.0-cudnn8-runtime-ubuntu18.04

RUN \
    apt-get update && \
    apt-get -y install wget \
    aria2 \
    rsync \
    cmake \
    build-essential \
    curl \
    git \
    vim \
    nvidia-cuda-toolkit

# Create environment
RUN conda env create -f /openfold/environment.yaml
SHELL ["conda", "run", "-n", "openfold-env", "/bin/bash", "-c"]
RUN pip install -r /openfold/requirements.txt
RUN pip install nvidia-pyindex \
    && pip install nvidia-dllogger
SHELL ["/bin/sh", "-c"]

#clone mmseqs2 and compile
RUN git clone https://github.com/soedinglab/MMseqs2.git
RUN mkdir -p /MMseqs2/build
WORKDIR /MMseqs2/build
RUN cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. ; make ; make install

# Download folding resources
RUN mkdir -p /openfold/resources
RUN wget -q -P /openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Certain tests need access to this file
RUN mkdir -p /openfold/tests/test_data/alphafold/common \
    && ln -rs /openfold/resources/stereo_chemical_props.txt /openfold/tests/test_data/alphafold/common

# Decompress test data
RUN gunzip /openfold/tests/test_data/sample_feats.pickle.gz

# Download pretrain alphafold weights
RUN /openfold/scripts/download_alphafold_params.sh /openfold/resources/

ENTRYPOINT ["/openfold/startup.sh"]
