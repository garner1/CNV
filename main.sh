#!/usr/bin/env bash

bam=$1				# full path to bam file
refgen=$2			# full path to reference genome
out=$3				# output dir
annotation=$4			# full path to annotation file refFlat.txt

rm -r $out/*
mkdir -p $out

cnvkit.py batch ~/Work/dataset/reduced_sequencing/{XZ68/outdata/GTCGTCGC.deduplicated.bam,XZ82/outdata/CATCATCC.deduplicated.bam,XZ83/outdata/TCACACGC.deduplicated.bam,XZ85/outdata/CTAACTCA.deduplicated.bam,XZ86/outdata/GAATCCGA.deduplicated.bam,XZ88/outdata/GTCGTTCC.deduplicated.bam,XZ89/outdata/CGTGTCGC.deduplicated.bam,XZ90/outdata/CGTGTGAG.deduplicated.bam} -n -m wgs -f $refgen --annotate $annotation --output-reference $out/my_flat_reference.cnn -d $out
