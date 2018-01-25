#!/usr/bin/env bash

bam=$1				# full path to bam file
baits=$2			# full path to baits bed file
refgen=$3			# full path to reference genome
out=$4				# output dir

mkdir -p $out
cnvkit.py batch $bam -n -t $baits -f $refgen --output-reference $out/my_flat_reference.cnn -d $out
