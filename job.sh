#!/bin/bash
#============ PBS Options ============
#QSUB -ug gr20001
#QSUB -q gr20001a
#QSUB -rn
#QSUB -W 1:00:00
#QSUB -A p=16:t=4:c=4:m=5G
#QSUB -e stderr.err
#QSUB -o stdout.out
#============ Shell Script ============

aprun -n $QSUB_PROCS -d $QSUB_THREADS -N $QSUB_PPN ./fdtd
