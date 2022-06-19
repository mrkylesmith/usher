#!/bin/bash
#
# Install 3seq 

git clone https://gitlab.com/mrkylesmith/3seq.git
cd 3seq
./3seq -g my3seqTable700 700
make
