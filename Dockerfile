ARG IMAGE_NAME=nvcr.io/nvidia/pytorch:21.06-py3

FROM ${IMAGE_NAME} AS base
ENV DEBIAN_FRONTEND=noninteractive


COPY . /workspace/openfold
WORKDIR /workspace/openfold

RUN pip install -r requirements_minimal.txt
RUN python setup.py install
