#!/bin/bash

# submitten mit:
# bash bci-markov-submit-wrapper-correct-full-long.sh 1000 analyse-1

ITERATIONS=$1
JOB_NAME_SUFFIX=$2

for grid in /data/idiv_chase/OliverPurschke/bci-input/grid.full.*.Rdata ; do
  for env in /data/idiv_chase/OliverPurschke/bci-input/env.*.Rdata ; do
    for columns in "1" "2" "3" "c(1,2)" "4" "c(4,5)" ; do

      case $env in 
	*.mean*)
          lenind=0
          ;;

	*.sd*)
          lenind=1
          ;;
      esac

      qsub -N bci.markov-$JOB_NAME_SUFFIX ~/bci/bci-markov-submit-correct-full-long.sh $grid $env $columns $lenind $ITERATIONS
    done
  done
done
