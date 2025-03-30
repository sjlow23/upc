#!/bin/bash

domain="$1"
use_assembly="$2"
species_taxid="$3"
outdir="$4"
metadir="$5"
offtargetdir="$6"
max_offtarget="$7"
database="$8"
assembly_level="$9"



# Download genomes using nt db
if [[ "$domain" == "viral" && "$use_assembly" == "no" ]]
then
	datasets download virus genome \
		taxon "$species_taxid" \
		--complete-only \
		--filename "$outdir"/offtarget.zip
	dataformat tsv virus-genome --force --package "$outdir"/offtarget.zip > "$metadir"/metadata_offtarget.tsv

	unzip "$outdir"/offtarget.zip -d "$offtargetdir"
	rm "$outdir"/offtarget.zip

	# If subsampling required
	if [[ "$max_offtarget" != "no" ]]
	then
		cat "$offtargetdir"/ncbi_dataset/data/genomic.fna | \
			seqkit sample -n "$subsample" -o "$outdir"/genomic.fna
			seqkit split --by-id --by-id-prefix "" -O "$offtargetdir" "$outdir"/genomic.fna
			rm "$outdir"/genomic.fna
			ls "$offtargetdir"/*.fna | awk -F "/" '{ print $NF }' > "$outdir"/offtarget_genomes_subsampled.txt
	else
		seqkit split --by-id --by-id-prefix "" -O "$offtargetdir" "$offtargetdir"/ncbi_dataset/data/genomic.fna
	fi

	rm -rf "$offtargetdir"/ncbi_dataset
	rm "$offtargetdir"/README.md "$offtargetdir"/md5sum.txt 

	# Generate offtarget genome list
	cat "$metadir"/metadata_offtarget.tsv | \
		awk -F "\t" '{ print $1 }' | \
		grep -v Accession | \
		awk -F "_" '{ print $1"_"$2".fna" }' | \
		sort | uniq > "$outdir"/offtarget_genomes.txt
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
		--filename "$outdir"/offtarget.zip
				
	unzip "$outdir"/offtarget.zip -d "$offtargetdir"

	dataformat tsv genome --package "$outdir"/offtarget.zip > "$metadir"/metadata_offtarget_assembly.tsv
	datasets summary virus genome taxon "$species_taxid" --as-json-lines | \
		dataformat tsv virus-genome \
		--fields accession,geo-location,geo-region,isolate-collection-date,host-common-name,host-name,virus-common-name,virus-name,virus-tax-id \
			> "$metadir"/metadata_offtarget.tsv
	
	cat "$metadir"/metadata_offtarget_assembly.tsv | \
		awk -F "\t" '{ print $1 }' | \
		grep -v Accession | \
		awk -F "_" '{ print $1"_"$2".fna" }' | \
		sort | uniq > "$outdir"/offtarget_genomes.txt


	# If subsampling required
	if [[ "$max_offtarget" != "no" ]]
	then	
		shuf -n "$max_offtarget" "$offtargetdir"/ncbi_dataset/fetch.txt > "$offtargetdir"/ncbi_dataset/tmp
		mv "$offtargetdir"/ncbi_dataset/tmp "$offtargetdir"/ncbi_dataset/fetch.txt
	fi

	datasets rehydrate --directory "$offtargetdir"
	rm "$outdir"/offtarget.zip

	mv "$offtargetdir"/ncbi_dataset/data/GC?_*/*.fna "$offtargetdir"/
	rm -rf "$offtargetdir"/ncbi_dataset "$offtargetdir"/README.md "$outdir"/assembly_accession.tmp

	find "$offtargetdir" -type f -name 'GC*.fna' -exec bash -c 'mv "$1" "$(dirname "$1")/$(basename "$1" | cut -d"_" -f1,2).fna"' _ {} \;

	if [[ "$max_offtarget" != "no" ]]
	then
		for genome in "$offtargetdir"/*.fna; do basename "$genome" >> "$outdir"/offtarget_genomes_subsampled.txt; done
	fi

fi	





