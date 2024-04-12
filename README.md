# wbp_scratch
Scratch scripts for performing analysis related to WormBase ParaSite data.

### `python3 parse_pdb.py pdb_dir`
**pdb_dir** should contain .pdb files of AlphaFold models. Each file is parsed using BioPython, and statistics derived from the extracted pLDDT scores. Prints CSV-formatted lines to the console of "pdb_filename,mean,median".