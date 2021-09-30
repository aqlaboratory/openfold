#!/bin/bash

FLAGS=""

while getopts ":v" option; do
    case $option in
        v)
            FLAGS=$(echo "-v $FLAGS" | xargs) # strip whitespace
            ;;
        *)
            echo "Invalid option: ${option}"
            ;;
    esac
done



python3 -m unittest $FLAGS "$@" || \
echo -e "\nTest(s) failed. Make sure you've installed all Python dependencies."
