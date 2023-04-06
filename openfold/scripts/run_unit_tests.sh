#!/bin/bash

CUDA_VISIBLE_DEVICES="0"

python3 -m unittest "$@" || \
echo -e "\nTest(s) failed. Make sure you've installed all Python dependencies."
