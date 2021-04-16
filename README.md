This repository hosts code (primarily R and Python scripts) and data related to the cold/heat stress response project.

Directory architecture:

`README.md`: (this file)
`src/`: data processing, statistical testing, visualization (R)
`nf.degA/`, `nf.degB/`, `nf.dmodA`, `nf.dmodB`, `nfb1-4`: nextflow pipeline configuration and results

Links to datasets:
- [Raw read count and normalized expression (Counts per Million) matrices for the inbred/hybrid experiment and the time-course inbred-only experiment]
- [List of DEGs at 1/25 hours after cold/heat treatment identified in the inbred/hybrid experiment]
- [List of DEGs showing variable response to stress among genotypes and their inferred cis/trans regulatory pattern]
- [Gene lists for each co-expression module identified from the time-course experiment]
- [Known and novel motifs discovered by motif mining for each set of DEGs under cold/heat stress as well as co-expression modules]
- [DEGs showing variable response to stress that are correctly predicted by the machine learning model (either the "B" or the "BMW_model")](https://s3.msi.umn.edu/zhoup-stress/71_share/08.variable.genes.tsv)

Links to R/Python scripts:
- [kmer.py](https://github.com/orionzhou/nf/blob/master/bin/kmer.py): kmer utilities, use `kmer.py -h` to find out more
  - `kmer.py locate`: find given kmers in sequence database and report locations
  - `kmer.py prepare_ml`: locate given kmers for given IDs in sequence db using various filters and prepare output for ML
  - `kmer.py getfasta`: extract fasta for given IDs in sequence db using various filters
- [fimo.py](https://github.com/orionzhou/nf/blob/master/bin/mmm/fimo.py): wrapper around `FIMO` from the [meme-suite](https://meme-suite.org/meme/)
- [streme.py](https://github.com/orionzhou/nf/blob/master/bin/mmm/streme.py): wrapper around `STRME` from the [meme-suite](https://meme-suite.org/meme/), output a meme file with found motifs and a tabular file with the exact motif locations
- [merge.fimo.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/merge.fimo.R): read multiple FIMO outputs and save as a tibble in R
- [merge.dreme.kmer.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/merge.dreme.kmer.R): read multiple meme outputs after running DREME/STRME and save as a tibble in R
- [merge.dreme.fimo.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/merge.dreme.fimo.R): read multiple motif location outputs after running DREME/STRME and save as a tibble in R
- [merge.dreme.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/merge.dreme.R): read multiple DREME outputs save as a tibble in R
- [ml_classification.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/ml_classification.R): train a machine learning model using RF/XGB/SVM algorithm, specifying holdout proportion, with down-sampling, using cross-validation, grid searching for hyperparameters in parallel, for detailed usage run `ml_classification.R --help`
- [ml_predict.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/ml_predict.R): read a trained model and predict outcome of a new dataset, for detailed usage run `ml_predict.R --help`


