#!/usr/bin/env bash

set -euo pipefail

FDS="$PWD/../build/fds"
TEST_DIR=.

clean_test() {
  rm -rf */output_*
}

run_test() {
  mkdir -p $TEST_DIR/$1/output_$2
  cd $TEST_DIR/$1/output_$2
  $FDS ../$2.fds

  # remove irrelevant files (they will differ since they contain timer info)
  rm *_cpu.csv
  rm *_steps.csv
  rm *.out

  # diff with the expected files
  if [ ! -z "$(diff -r . ../$2_expected)" ]; then
    echo "ERROR: $1/$2"
  fi
}

run_config_file() {
  dir=$(echo $1 | awk -F'/' '{ print $1 }')
  file=$(echo $1 | awk -F'/' '{ print $2 }' | sed 's/.fds//g')

  run_test $dir $file
}

if [ $# -gt 0 ]; then
  case $1 in
    clean) clean_test;;
    all) echo "TODO: run all tests" ;;
    *)
      run_config_file $1
      ;;
  esac
else
  # TODO: create a help message
  echo "Error: expected parameters"
fi
