rule download_hg:
	input:
		#infile = "myseqs.fasta"
	
	output:
		status = OUTDIR + "status/download_hg.txt"

	params:
		db = DATABASE,

	threads: 4

	shell:
		"""
		if [ ! -s {params.db}/GCF_000001405.40_GRCh38.p14_genomic.fna ] 
		then
			wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz -P {params.db}
			gunzip -f {params.db}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
		fi

		touch {output.status}
		"""
