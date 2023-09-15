#!/bin/bash
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

echo "Setting up Python3 from CVMFS..."
# setup python v3_7_2
# setup root v6_18_04 -q e17:prof:py3

setup python v3_9_13
setup root v6_26_06b -q e20:p3913:prof

if [ $? -ne 0 ]; then
    echo "Setup failed. Script must be sourced, not run."
    echo "Also check if the CVMFS environment is correct."
fi

source ./env/bin/activate
