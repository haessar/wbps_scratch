# wbp_scratch
Scratch scripts for performing analysis related to WormBase ParaSite data.

### `./parse_pdb.py pdb_dir`
**pdb_dir** should contain .pdb files of AlphaFold (AF) models. Each file is parsed using BioPython, and statistics derived from the extracted pLDDT scores. Prints CSV-formatted lines to the console of "pdb_filename,mean,median,stdev,var,max,min,perc_confident", where *perc_confident* is defined as the percentage of AF model residues that are above the 70% ("Confident") threshold.

### `pdb_analysis.ipynb`
scp (large) datasets from `MARS:/users/whh2g/sharedscratch/parse_pdb/` and place in `data/from_MARS/`. Uses `matplotlib` to produce plots in output dir `plots/`.

### `./get_all_species.py`
Download all gzipped gff3 and genomic fasta files for all species from the current WBP release, and store in `data/from_WBP/`.
