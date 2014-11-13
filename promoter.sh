#!/bin/bash
#PBS -N promoter_prediction
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=12
#PBS -l walltime=32:00:00
#PBS -d /home/mrals/promoter
#------------------------------------------------
# Title: promoter.sh

CORES=12
export TARGET=CAC.txt
export PVAL=0.01
# Minimum and maximum distance between bipartite motifs
MINDIST=10
MAXSIT=30
# DIRECTORIES
INDIR=motifs
PSPM=pspm
export OUTDIR=rawresults
RESULTS=results
#./promoter_prep.rb INDIR PSPM


# Process the promoter prediction in parallel

#parallel -j$CORES 'mast {} $TARGET -o $OUTDIR/{/.} -hit_list -mt $PVAL -comp > $OUTDIR/{/.}.txt' ::: `/usr/bin/ls $PSPM/*`


# Process the results

./promoter_parse.rb $OUTDIR $RESULTS $CORES $MINDIST $MAXDIST
