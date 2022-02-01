FROM nvidia/cuda:11.0-cudnn8-runtime-ubuntu18.04

# I'm not sure why i needed both opencl and cuda here, but the relax phase of the script needed opencl
RUN apt-get update && apt-get install -y wget cuda-minimal-build-11-0 nvidia-opencl-dev
RUN wget -P /tmp \
    "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH

COPY environment.yml /opt/openfold/environment.yml

# this needs to be run separately so that nvidia-dllogger will install properly
RUN /opt/conda/bin/python3 -m pip install nvidia-pyindex

# installing into the base environment since the docker container wont do anything other than run openfold
RUN conda env update -n base --file /opt/openfold/environment.yml && conda clean --all

COPY openfold /opt/openfold/openfold
COPY scripts /opt/openfold/scripts
COPY run_pretrained_openfold.py /opt/openfold/run_pretrained_openfold.py
COPY train_openfold.py /opt/openfold/train_openfold.py
COPY setup.py /opt/openfold/setup.py
COPY lib/openmm.patch /opt/openfold/lib/openmm.patch
RUN wget -q -P /opt/openfold/openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
RUN patch -p0 -d /opt/conda/lib/python3.7/site-packages/ < /opt/openfold/lib/openmm.patch
RUN python3 /opt/openfold/setup.py install
