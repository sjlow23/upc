#!/bin/bash

querylist=$1
targetdir=$2
alndir=$3
segmentdir=$4
grouped=$5
genomestatus=$6
primerexpand=$7
mode=$8
thread=$9


if [[ "$mode" == "primer" ]]; then
	while read primer
	do
		touch "$alndir"/keep_"$primer".txt

		# Extract sequences with "not amplified" and filter by primer, then save to keep_"$primer".txt
		grep "Not amplified" $genomestatus | grep -w "$primer" | cut -f1 > "$alndir"/keep_"$primer".txt

		# If file is not empty
		if [[ -s "$alndir"/keep_"$primer".txt ]]; then

			# Get perfect genomes from same primer_set
			grep -w "$primer" $grouped | \
				awk -F "\t" '$4 == 0 && $5 == 0' | \
				sort -k2,3 | \
				head -n1 | \
				cut -f1,3 > "$alndir"/keep_"$primer"_perfect.txt
			
			primerset=$(cut -f2 "$alndir"/keep_"$primer"_perfect.txt | head -n1)

			# Get fwd and rev primer
			fwd_primer=$(grep -w "$primerset" $primerexpand | cut -f3)
			rev_primer=$(grep -w "$primerset" $primerexpand | cut -f4)

			for i in `cut -f1 "$alndir"/keep_"$primer"_perfect.txt`; do 
				cat "$targetdir"/"$i".fna | awk '{ print $1 }' | sed 's/>/>REF_/g' >> "$alndir"/check_"$primer".fasta
			done

			## For segmented genomes, need to ensure only contig of interest is extracted
			## Need to know exactly which segment is targeted, and only extract those segments for seq comparison
			for i in `cat "$alndir"/keep_"$primer".txt`; do 
				cat "$targetdir"/"$i".fna | awk '{ print $1 }' >> "$alndir"/check_"$primer".fasta
			done

			# Align the sequences using mafft
			mafft --thread $thread "$alndir"/check_"$primer".fasta > "$alndir"/check_"$primer".aln

			# Locate primer sequences in alignment
			fwd_coord=$(seqkit seq --upper-case "$alndir"/check_"$primer".aln | seqkit locate -i -p "$fwd_primer" | cut -f5-6 | grep -v start | sort | uniq | head -1 | sed 's/\t/:/g') 
			rev_coord=$(seqkit seq --upper-case "$alndir"/check_"$primer".aln | seqkit locate -i -p "$rev_primer" | cut -f5-6 | grep -v start | sort | uniq | head -1 | sed 's/\t/:/g')

			cat "$alndir"/check_"$primer".aln | seqkit subseq -r $fwd_coord > "$alndir"/check_"$primer"_fwd.fasta
			cat "$alndir"/check_"$primer".aln | seqkit subseq -r $rev_coord | seqkit seq --reverse --complement > "$alndir"/check_"$primer"_rev.fasta

			# Plot msa plot
			Rscript scripts/plot_msaprimer.R \
				"$alndir"/check_"$primer"_fwd.fasta \
				"$alndir"/check_"$primer"_rev.fasta \
				"$alndir"/"$primer"_fwd.png \
				"$alndir"/"$primer"_rev.png \
				"$primer"

			magick convert "$alndir"/"$primer"_fwd.png  -bordercolor white -border 250x0 "$alndir"/"$primer"_rev.png +append "$alndir"/"$primer".png

			# Clean up temporary files
			#rm "$alndir"/keep_"$primer".txt "$alndir"/keep_"$primer"_perfect.txt 
			#rm "$alndir"/"$primer"_fwd.png "$alndir"/"$primer"_rev.png
			#rm "$alndir"/check_"$primer"_fwd.fasta "$alndir"/check_"$primer"_rev.fasta

		else
			continue
		fi
	done < $querylist
fi


if [[ "$mode" == "probe" ]]; then
	while read probe
	do
		touch "$alndir"/probe_"$probe".txt

		# Extract sequences with "not amplified" and filter by probe, then save 
		grep "Not present" $genomestatus | grep -w "$probe" | cut -f1 > "$alndir"/probe_"$probe".txt

		# If file is not empty
		if [[ -s "$alndir"/probe_"$probe".txt ]]; then

			# Get perfect genomes from same probe set
			grep -w "$probe" $grouped | sort -nk4 -k2 | head -1 | cut -f1-2 > "$alndir"/probe_"$probe"_perfect.txt
			
			probeset=$(cut -f2 "$alndir"/probe_"$probe"_perfect.txt | head -n1)

			# Get probe seq
			probeseq=$(grep -w "$probeset" $primerexpand | cut -f3)

			for i in `cut -f1 "$alndir"/probe_"$probe"_perfect.txt`; do 
				cat "$targetdir"/"$i".fna | awk '{ print $1 }' | sed 's/>/>REF_/g' >> "$alndir"/probe_"$probe".fasta
			done

			for i in `cat "$alndir"/probe_"$probe".txt`; do 
				cat "$targetdir"/"$i".fna | awk '{ print $1 }' >> "$alndir"/probe_"$probe".fasta
			done

			# Align the sequences using mafft
			mafft --thread $thread "$alndir"/probe_"$probe".fasta > "$alndir"/probe_"$probe".aln

			# Locate probe in alignment
			probe_coord=$(seqkit seq --upper-case "$alndir"/probe_"$probe".aln | seqkit locate -i -p "$probeseq" | cut -f5-6 | grep -v start | sort | uniq | head -1 | sed 's/\t/:/g') 
			
			cat "$alndir"/probe_"$probe".aln | seqkit subseq -r $probe_coord > "$alndir"/probe_"$probe".fasta
			
			# Plot msa plot
			Rscript scripts/plot_msaprobe.R \
				"$alndir"/probe_"$probe".fasta \
				"$alndir"/"$probe".png \
				"$probe"

			# Clean up temporary files
			#rm "$alndir"/probe_"$probe".txt "$alndir"/probe_"$probe"_perfect.txt 

		fi
	done < $querylist
fi
