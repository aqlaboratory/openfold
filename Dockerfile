ARG IMAGE_NAME=nvcr.io/nvidia/pytorch:21.06-py3

FROM ${IMAGE_NAME} AS base
ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && \
    apt install -y aria2

COPY . /workspace/openfold
WORKDIR /workspace/openfold

RUN pip install -r requirements_minimal.txt
RUN python setup.py install

# TODO add all dependencies needed for inference

RUN wget -q -P openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt && \
    mkdir -p tests/test_data/alphafold/common && \
    ln -rs openfold/resources/stereo_chemical_props.txt tests/test_data/alphafold/common

RUN scripts/download_alphafold_params.sh openfold/resources

RUN gunzip tests/test_data/sample_feats.pickle.gz
