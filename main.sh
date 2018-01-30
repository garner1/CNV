#!/usr/bin/env bash

bam=$1				# full path to bam file
refgen=$2			# full path to reference genome
out=$3				# output dir
annotation=$4			# full path to annotation file refFlat.txt

mkdir -p $out

cnvkit.py batch $bam -n -m wgs -f $refgen --annotate $annotation --output-reference $out/my_flat_reference.cnn -d $out
