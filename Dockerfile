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

# install miniconda
RUN \
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh && \
    sh miniconda.sh -b -p /tools/miniconda && \
    rm miniconda.sh && \
    /tools/miniconda/bin/conda init && \
    /tools/miniconda/bin/conda clean -a
ENV PATH "$PATH:/tools/miniconda/bin"

COPY . /openfold
WORKDIR /openfold

# Create environment
RUN conda env create -f environment.yaml
SHELL ["conda", "run", "-n", "openfold-env", "/bin/bash", "-c"]
RUN pip install -r requirements.txt
RUN pip install nvidia-pyindex \
    && pip install nvidia-dllogger
SHELL ["/bin/sh", "-c"]

#clone mmseqs2 and compile
WORKDIR /
RUN git clone https://github.com/soedinglab/MMseqs2.git
RUN mkdir -p /MMseqs2/build
WORKDIR /MMseqs2/build
RUN cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. ; make ; make install
WORKDIR /openfold

# Download folding resources
RUN mkdir -p resources
RUN wget -q -P resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Download pretrain alphafold weights
RUN scripts/download_alphafold_params.sh /openfold/resources/

# execute startup script
chmod +x /openfold/startup.sh
ENTRYPOINT ["startup.sh"]
