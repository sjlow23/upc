#!/bin/bash

#$1: genome
#$2: target genome dir
#$3: bbmap_amplicons dir
#$4: target_nohits
#$5: primer name

msa.sh in="$2"/"$1".fna \
	noheader=t \
	trimreaddescriptions=t \
	out="$3"/"$1"_"$5".fwd.sam \
	ref="$4"/"$5".fwd.fasta cutoff=0.75

msa.sh in="$2"/"$1".fna \
	noheader=t \
	trimreaddescriptions=t \
	addr=t \
	out="$3"/"$1"_"$5".rev.sam \
	ref="$4"/"$5".rev.fasta cutoff=0.75

cutprimers.sh \
	trimreaddescriptions=t \
	sam1="$3"/"$1"_"$5".fwd.sam \
	sam2="$3"/"$1"_"$5".rev.sam \
	in="$2"/"$1".fna \
	out="$3"/"$1"--"$5".fasta \
	include=t \
	fake=f

filename=$(basename "$3"/"$1"--"$5".fasta .fasta)
awk -v name="$filename" '/^>/ {print ">"name} !/^>/ {print}' "$3"/"$1"--"$5".fasta > "$3"/"$1"--"$5".fasta.2
mv "$3"/"$1"--"$5".fasta.2 "$3"/"$1"--"$5".fasta

if [[ -s "$3"/"$1"--"$5".fasta ]]; then
	for i in "$3"/"$1"_"$5".fwd.sam "$3"/"$1"_"$5".rev.sam; do \
		awk -F "\t" 'BEGIN {
						comp["A"] = "T"; comp["T"] = "A"; comp["C"] = "G"; comp["G"] = "C"
					}
					$2 != 4 && $2 != 20 {
						if ($2 ~ /16/) {
							# If the alignment is on the reverse strand (FLAG 16), reverse complement the sequence
							seq = $10
							rev_comp = ""
							for (i = length(seq); i > 0; i--) {
								rev_comp = rev_comp comp[substr(seq, i, 1)]
							}
							$10 = rev_comp
						}
						# Print the line with the sequence either reversed or not
						print
					}' OFS="\t" "$i" | sed 's/^r_//g' > "$i".rc
	done

	awk -F "\t" '{ print $3"--"$1, $1, $4, $10 }' OFS="\t" "$3"/"$1"_"$5".fwd.sam.rc | sort -k1,1 > "$3"/"$1"_"$5".fwd
	awk -F "\t" '{ print $3"--"$1, $1, $4+length($10)-1, $10 }' OFS="\t" "$3"/"$1"_"$5".rev.sam.rc | sort -k1,1 >> "$3"/"$1"_"$5".rev


	join -1 1 -2 1 -a 1 -a 2 -e NA -o 1.1,1.2,2.2,1.3,2.3,1.4,2.4 "$3"/"$1"_"$5".fwd "$3"/"$1"_"$5".rev | \
			sed 's/--/\t/1' | \
			awk '{ print $1"--"$2, $1":"$5"+"$6, $2, "X", $7, $8 }' OFS="\t" | \
			sed 's/\t/ /g' | \
			sed 's/ /\t/1' > "$3"/"$1"_"$5".header

	seqkit replace --kv-file "$3"/"$1"_"$5".header --pattern "^(\\S+)" --replacement "{kv}" --out-file "$3"/"$1"_"$5"_amp.fasta "$3"/"$1"--"$5".fasta
	rm "$3"/"$1"_"$5".fwd.sam.rc "$3"/"$1"_"$5".rev.sam.rc "$3"/"$1"_"$5".fwd "$3"/"$1"_"$5".rev "$3"/"$1"_"$5".header "$3"/"$1"_"$5".fwd.sam "$3"/"$1"_"$5".rev.sam "$3"/"$1"--"$5".fasta
else
	rm "$3"/"$1"--"$5".fasta "$3"/"$1"_"$5".fwd.sam "$3"/"$1"_"$5".rev.sam
fi



