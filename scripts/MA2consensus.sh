#!/bin/bash

mkdir ./analysis/separated_strains
mkdir ./analysis/consensus
echo "decompressing input ..."
gunzip ./analysis/MAfile.fa.gz

echo "counting amino acid frequencies ..."
python3 aaFreqs.py

echo "making consensuses ..."
Rscript consensus.R

echo "Process completed."
