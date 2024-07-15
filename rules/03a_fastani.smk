rule fastani_target_genomes:
	input:
		inputdir = GENOMES_TARGET,
		genomelist = OUTDIR + "genomes_target.txt"
	output:
		fastani = OUTDIR + "phylogeny/fastani.tsv",
		dist = OUTDIR + "phylogeny/distance.tsv",
		treefile = OUTDIR + "phylogeny/combined_target.nwk",
		status = OUTDIR + "status/phylogeny.txt"
	conda: "../envs/fastani.yaml"
	params:
		outdir = OUTDIR + "phylogeny",
	threads: 16
	shell:
		"""
		genomedir="{input.inputdir}"
		dirname=$(echo $genomedir | sed 's/\/\$//g')
		sed "s#^#$dirname\/#g" {input.genomelist} | sed "s/$/.fna/g" > {params.outdir}/genomes.txt
		fastANI --fragLen 200 --queryList {params.outdir}/genomes.txt --refList {params.outdir}/genomes.txt -t {threads} -o {output.fastani}

		sed -i "s#$dirname\/##g" {output.fastani}
		awk -F "\\t" '{{ print $0, 1-($3*($4/$5/100)) }}' OFS="\\t" {output.fastani} | sed 's/.fna//g' | cut -f1-2,6 > {output.dist}
		
		Rscript ./scripts/phylo.R {output.dist} {output.treefile}
		rm {params.outdir}/genomes.txt
		touch {output.status}
		"""


rule tree_metadata:
	input:
		infile = rules.fastani_target_genomes.output.status,
	output:
		subset_meta = OUTDIR + "phylogeny/metadata_subset.tsv"
	params: 
		outdir = OUTDIR,
		genomedir = GENOMES_TARGET,
		metadata = OUTDIR + "metadata_target.tsv"
	shell:
		"""
		grep ">" {params.genomedir}/*.fna | awk '{{ print $1 }}' | \
			awk -F "/" '{{ print $NF }}' | \
			sed 's/.fna//g' | \
			sed 's/:/\\t/1' | \
			sed 's/>//g' | \
			awk '{{ print $2, $1 }}' OFS="\\t" > {params.outdir}/assembly_accession.txt 

		cut -f1,61,77,89,164 {params.metadata} > {output.subset_meta}
		"""
