#!/bin/bash

configsDir=$1
refDir=$2

# For configs dir containing files of form "params_X.config", each producing a reference dir "ref_X"
for x in $configsDir/*
do
    nextflow run reference_update.nf -c $x -with-trace -qs 8
    echo $x | sed -e "s/^configs\/params/ref/" -e "s/.config$//" | xargs mv -t $refDir/backup
done
