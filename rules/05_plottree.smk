rule plot_tree:
	input:
		"output_tree/denv2.nwk",
		"data/sequences_metadata.csv"
	output:
		"plots/denv2_linear.pdf",
		"plots/denv2_circular.pdf",
		"plots/denv2_facet.pdf"
	script:
		"scripts/plot_tree.R"