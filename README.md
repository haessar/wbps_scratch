# wbp_scratch
Scratch scripts for performing analysis related to WormBase ParaSite data.

### `./parse_pdb.py pdb_dir`
**pdb_dir** should contain .pdb files of AlphaFold (AF) models. Each file is parsed using BioPython, and statistics derived from the extracted pLDDT scores. Prints CSV-formatted lines to the console of "pdb_filename,mean,median,stdev,var,max,min,perc_confident", where *perc_confident* is defined as the percentage of AF model residues that are above the 70% ("Confident") threshold.

### `pdb_analysis.ipynb`
scp (large) datasets from `MARS:/users/whh2g/sharedscratch/parse_pdb/` and place in `data/from_MARS/`. Uses `matplotlib` to produce plots in output dir `plots/`.

### `./get_all_species.py output_dir`
Download all gzipped gff3 and genomic fasta files for all species from the current WBP release, and store in **output_dir**.

### `./prepare_geenuff_inputs.py input_dir output_dir`
Create necessary directory structure to be compatable with the import2geenuff.py script from GeenuFF. **input_dir** should contain the raw gzipped gff3 and genomic fasta files as gathered by `get_all_species.py`.

### `determine_training_species.ipynb`
Uses TSV of BUSCO scores for annotation and assembly of all WBP species to select candidates for reannotation, model training and validation in Helixer. Outputs text lists of candidate species to `data/helixer_training_species/`.

### `./prepare_helixer_training_inputs.py train_set_file valid_set_file train_dir`
Prepare symbolic links in **train_dir** with the following hierarchy, using **train_set_file** and **valid_set_file** species lists:
```
<train_dir>
├── training_data.species_01.h5 -> ../h5s/speciesA/test_data.h5
├── training_data.species_02.h5 -> ../h5s/speciesB/test_data.h5
├── validation_data.species_03.h5 -> ../h5s/speciesC/test_data.h5
└── validation_data.species_04.h5 -> ../h5s/speciesD/test_data.h5
```

### `./get_gene_map.py gffcmp_refmap`
Extract gene IDs from given GFFCompare .refmap file **gffcmp_refmap** using a regex pattern. Prints TSV-formatted lines to console of `ref_gene_id  target_gene_id`.

### `./busco_preprocessing.sh seqfile annfile`
Sort annotation GFF3 **annfile** and translate CDS to protein sequences using genome fasta **seqfile**. These can then be used for running BUSCO in "proteins" mode.

### `./get_prot_seq_for_uniprot_acc.py uniprot_acc`
Call the EBI AlphaFold API to get the protein sequence for a given **uniprot_acc**.

### `./find_microexon_genes.py output_dir`
Prints to console genes that comply with criteria defined within script, followed by the total number. If (optional) **output_dir** is specified, write MD-formatted exon-lengths for each microexon gene.
