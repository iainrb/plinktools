#! /usr/bin/env bash

# workaround for annoying NFS error in deleting directories within Python tests

PYTHON_TEST='./test.py'

if [ ! -f $PYTHON_TEST ]; then
    echo "$PYTHON_TEST does not exist; should run from the plinktools source code directory"
    exit 1
else
    python $PYTHON_TEST
    echo "Finished running $PYTHON_TEST";
    rm -Rf output_test_*
fi