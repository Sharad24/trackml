#!/bin/bash
#
# This script tests whether the output of a certain ACTS framework example is
# reproducible between single-threaded and multi-threaded runs. For example,
# "./testReproducibility.sh Extrapolation" will run the ExtrapolationExample
# in single-threaded and multi-threaded mode and check whether the output is
# the same aside from threading-induced event reordering.
#
set -uo pipefail

# Check whether the user did specify the name of the example to be run
ARGC=$#
if [[ $ARGC -lt 2 ]]; then
  echo ""
  echo " Usage: "$0" <example> <nevents> <output1> [<output2> ...]"
  echo ""
  echo " <example> is the example name (which is the executable name without the leading 'ACTFW' and the trailing 'Example')"
  echo " <outputN> is the output name (which is the output file name without the trailing '.root')"
  echo ""
  exit 42
fi

# Compute the name of the example executable
executable="ACTFW$1Example -n $2 --output-root true"

# Compute the output file names
for ((i = 3; i <= $ARGC; i++)); do
  eval output=\$${i}.root
  eval outputs[$i]=$output
done

# Drop any remaining output file from previous runs of the example
for output in "${outputs[@]}"; do
  rm -f $output ST$output MT$output
done

# Run the example in multi-threaded mode
eval "${executable}"
result=$?
if [[ result -ne 0 ]]; then
  echo "Multi-threaded run failed!"
  exit $result
fi

# Back up the multi-threaded results
for output in "${outputs[@]}"; do
  mv $output MT$output
done

# Run the example in single-threaded mode
eval "ACTSFW_NUM_THREADS=1 ${executable}"
result=$?
if [[ result -ne 0 ]]; then
  echo "Single-threaded run failed!"
  exit $result
fi

# Back up the single-threaded results
for output in "${outputs[@]}"; do
  mv $output ST$output
done

# Check whether the results were identical (up to thread-induced reordering)
for output in "${outputs[@]}"; do
  # Compare the active results files
  cmd="root -b -q -l -x -e '.x compareRootFiles.cpp(\"ST$output\", \"MT$output\")'"
  eval $cmd
  result=$?

  # If the results were different, abort and return the output status code
  if [[ result -ne 0 ]]; then
    exit $result
  fi

  # Otherwise clean up and continue
  rm ST$output MT$output
done
