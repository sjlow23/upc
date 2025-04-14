#!/bin/bash

domain="$1"
use_assembly="$2"
species_taxid="$3"
outdir="$4"
metadir="$5"
targetdir="$6"
max_target="$7"
min_genome_size="$8"
database="$9"
assembly_level="${10}"



# Download genomes using nt db
if [[ "$domain" == "viral" && "$use_assembly" == "no" ]]
then
	datasets download virus genome \
		taxon "$species_taxid" \
		--complete-only \
		--filename "$outdir"/target.zip
	dataformat tsv virus-genome --force --package "$outdir"/target.zip > "$metadir"/metadata_target.tsv

	unzip "$outdir"/target.zip -d "$targetdir"
	rm "$outdir"/target.zip

	# If subsampling required
	if [[ "$max_target" != "no" ]]
	then
		cat "$targetdir"/ncbi_dataset/data/genomic.fna | \
			seqkit seq --min-len "$min_genome_size" | \
			seqkit sample -n "$max_target" -o "$targetdir"/genomic.fna
			seqkit split --by-id --by-id-prefix "" -O "$targetdir" "$targetdir"/genomic.fna
			rm "$targetdir"/genomic.fna
			ls "$targetdir"/*.fna | awk -F "/" '{ print $NF }' > "$outdir"/target_genomes_subsampled.txt
	else
		seqkit split --by-id --by-id-prefix "" -O "$targetdir" "$targetdir"/ncbi_dataset/data/genomic.fna
	fi

	rm -rf "$targetdir"/ncbi_dataset
	rm "$targetdir"/README.md "$targetdir"/md5sum.txt 

	# Remove genomes smaller than min genome size
	for genome in $targetdir/*.fna; do
		if [[ `seqkit stats $genome | tail -1 | awk '{ print $7 }' | sed 's/,//g'` -lt "$min_genome_size" ]]
		then
			rm $genome
		fi
	done

	# Generate target genome list
	awk -F "\t" '$27 >= '"$min_genome_size"'' "$metadir"/metadata_target.tsv | \
		cut -f1 | \
		grep -v Accession | \
		awk -F "\t" '{ print $1".fna" }' | \
		sort | uniq > "$outdir"/target_genomes.txt
fi


# Download genomes using assemblies db
if [[ "$domain" == "viral" && "$use_assembly" == "yes" ]]
then
	datasets download genome taxon "$species_taxid" \
		--assembly-source "$database" \
		--assembly-version latest \
		--assembly-level "$assembly_level" \
		--include genome,seq-report \
		--dehydrated \
		--filename "$outdir"/target.zip
				
	unzip "$outdir"/target.zip -d "$targetdir"

	dataformat tsv genome --package "$outdir"/target.zip > "$metadir"/metadata_target_assembly.tsv
	datasets summary virus genome taxon "$species_taxid" --as-json-lines | \
		dataformat tsv virus-genome \
		--fields accession,geo-location,geo-region,isolate-collection-date,host-common-name,host-name,virus-common-name,virus-name,virus-tax-id \
			> "$metadir"/metadata_target.tsv
	
	awk -F "\t" '$144 >= "'$min_genome_size'"' "$metadir"/metadata_target_assembly.tsv | \
		cut -f1 | \
		grep -v Accession | \
		awk -F "_" '{ print $1"_"$2".fna" }' | \
		sort | uniq > "$outdir"/target_genomes.txt


	# If subsampling required
	if [[ "$max_target" != "no" ]]
	then	
		shuf -n "$max_target" "$targetdir"/ncbi_dataset/fetch.txt > "$targetdir"/ncbi_dataset/tmp
		mv "$targetdir"/ncbi_dataset/tmp "$targetdir"/ncbi_dataset/fetch.txt
	fi

	datasets rehydrate --directory "$targetdir"
	rm "$outdir"/target.zip
		
	# for i in "$targetdir"/ncbi_dataset/data/GC?_*/sequence_report.jsonl; do dataformat tsv genome-seq --inputfile $i | \
	# 	awk -F "\\t" '{ print $7, $1 }' OFS="\\t" | \
	# 	grep -v Accession >> "$outdir"/assembly_accession.tmp; done
	# sort -k1b,1 "$outdir"/assembly_accession.tmp > "$outdir"/assembly_accession.txt

	mv "$targetdir"/ncbi_dataset/data/GC?_*/*.fna "$targetdir"/
	rm -rf "$targetdir"/ncbi_dataset "$targetdir"/README.md #"$outdir"/assembly_accession.tmp

	find "$targetdir" -type f -name 'GC*.fna' -exec bash -c 'mv "$1" "$(dirname "$1")/$(basename "$1" | cut -d"_" -f1,2).fna"' _ {} \;

	for genome in "$targetdir"/*.fna; do
		if [[ `seqkit stats $genome | tail -1 | awk '{ print $7 }' | sed 's/,//g'` -lt "$min_genome_size" ]]
		then
			rm $genome
		fi
	done

	if [[ "$max_target" != "no" ]]
	then
		for genome in "$targetdir"/*.fna; do basename "$genome" >> "$outdir"/target_genomes_subsampled.txt; done
	fi

	# Make lookup file for assembly and nt accessions
	for genome in "$targetdir"/*.fna; do
		while read -r line; do
			header=$(echo "$line" | awk '{print $1}' | sed 's/^>//')
			echo -e "$header\t$(basename "$genome")" >> "$outdir"/target_assembly_accession.txt.tmp
		done < <(grep "^>" "$genome")
	done
	sort -k1,1 "$outdir"/target_assembly_accession.txt.tmp | sed 's/.fna//g' > "$outdir"/target_assembly_accession.txt
	rm "$outdir"/target_assembly_accession.txt.tmp

fi	





