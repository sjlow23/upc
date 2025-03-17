#!/bin/bash


msa.sh in="$2"/"$1".fna \
	noheader=t \
	trimreaddescriptions=t \
	out="$3"/"$1".fwd.sam \
	ref="$4"/"$5".fwd.fasta cutoff=0.75

msa.sh in="$2"/"$1".fna \
	noheader=t \
	trimreaddescriptions=t \
	addr=t \
	out="$3"/"$1".rev.sam \
	ref="$4"/"$5".rev.fasta cutoff=0.75

cutprimers.sh \
	trimreaddescriptions=t \
	sam1="$3"/"$1".fwd.sam \
	sam2="$3"/"$1".rev.sam \
	in="$2"/"$1".fna \
	out="$3"/"$1".fasta \
	include=t \
	fake=f

#rm "$3"/"$1".fwd.sam* "$3"/"$1".rev.sam