## This repository hosts code (primarily R and Python scripts) and data related to the cold/heat stress response project.

### Directory structure:
- `README.md`: (this file)
- `src/`: data processing, statistical testing, visualization (R)
- `nf.degA/`, `nf.degB/`, `nf.dmodA`, `nf.dmodB`, `nfb1-4`: nextflow pipeline configuration and results

### Links to datasets:
- [Raw read count and normalized expression (Counts per Million) matrices for the inbred/hybrid experiment and the time-course inbred-only experiment]()
- [List of DEGs at 1/25 hours after cold/heat treatment identified in the inbred/hybrid experiment]()
- [List of DEGs showing variable response to stress among genotypes and their inferred cis/trans regulatory pattern]()
- [Gene lists for each co-expression module identified from the time-course experiment]()
- [Known and novel motifs discovered by motif mining for each set of DEGs under cold/heat stress as well as co-expression modules]()
- [DEGs showing variable response to stress that are correctly predicted by the machine learning model (either the "B" or the "BMW_model")](https://s3.msi.umn.edu/zhoup-stress/71_share/08.variable.genes.tsv)

### Links to R/Python scripts:
- [kmer.py](https://github.com/orionzhou/nf/blob/master/bin/kmer.py): kmer utilities, use `kmer.py -h` to find out more
'''
    usage: kmer.py [-h] {locate,prepare_ml,getfasta} ...

    kmer utilities

    optional arguments:
      -h, --help            show this help message and exit

    available commands:
      {locate,prepare_ml,getfasta}
        locate              find given kmers in sequence database and report
                            locations
        prepare_ml          locate given kmers for given IDs in sequence db using
                            various filters and prepare output for ML
        getfasta            extract fasta for given IDs in sequence db using
                            various filters
'''
- [fimo.py](https://github.com/orionzhou/nf/blob/master/bin/mmm/fimo.py): wrapper around `FIMO` from the [meme-suite](https://meme-suite.org/meme/)

    usage: fimo.py [-h] {locate,filter,bed2wide,prepare_ml} ...
    
    fimo utilities
    
    optional arguments:
      -h, --help            show this help message and exit
    
    available commands:
      {locate,filter,bed2wide,prepare_ml}
        locate              run fimo to find given motifs in input sequences
        filter              filter BED file using window size / epigenetic marks
        bed2wide            convert BED file to machine learing tables
        prepare_ml          pipeline to find motifs and output in BED / ML input
                            table

- [streme.py](https://github.com/orionzhou/nf/blob/master/bin/mmm/streme.py): wrapper around `STRME` from the [meme-suite](https://meme-suite.org/meme/), output a meme file with found motifs and a tabular file with the exact motif locations

    usage: streme.py [-h] {addscore,xml2tsv,pipe} ...
    
    STREME utilities
    
    optional arguments:
      -h, --help            show this help message and exit
    
    available commands:
      {addscore,xml2tsv,pipe}
        addscore            add score_thresh to STREME output
        xml2tsv             convert STREME xml output to tsv
        pipe                run STREME pipeline

- [merge.fimo.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/merge.fimo.R): read multiple FIMO outputs and save as a tibble in R
- [merge.dreme.kmer.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/merge.dreme.kmer.R): read multiple meme outputs after running DREME/STRME and save as a tibble in R
- [merge.dreme.fimo.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/merge.dreme.fimo.R): read multiple motif location outputs after running DREME/STRME and save as a tibble in R
- [merge.dreme.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/merge.dreme.R): read multiple DREME outputs save as a tibble in R
- [ml_classification.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/ml_classification.R): train a machine learning model using RF/XGB/SVM algorithm, specifying holdout proportion, with down-sampling, using cross-validation, grid searching for hyperparameters in parallel, for detailed usage run `ml_classification.R -h`

    usage: /home/springer/zhoux379/git/nf/bin/mmm/ml_classification.R
           [-h] [--perm PERM] [--alg ALG] [--holdout HOLDOUT] [--fold FOLD]
           [--fold_repeat FOLD_REPEAT] [--nlevel NLEVEL] [--downsample]
           [--seed SEED] [--response RESPONSE] [--cpu CPU]
           fi fo1
    
    Run machine learning classification on given dataset
    
    positional arguments:
      fi                    input dataset
      fo1                   output metrics file
    
    optional arguments:
      -h, --help            show this help message and exit
      --perm PERM           number permutations [default: 1]
      --alg ALG             ML algorithm [default: rf]
      --holdout HOLDOUT     proportion data to hold out for test [default: 0.8]
      --fold FOLD           cv fold [default: 5]
      --fold_repeat FOLD_REPEAT
                            repeat in each cv [default: 1]
      --nlevel NLEVEL       levels of hyperparameters to tune [default: 3]
      --downsample          downsample to get balanced [default: False]
      --seed SEED           random seed [default: 26]
      --response RESPONSE   response variable name [default: status]
      --cpu CPU             num. processors to use [default: 1]

- [ml_predict.R](https://github.com/orionzhou/nf/blob/master/bin/mmm/ml_predict.R): read a trained model and predict outcome of a new dataset, for detailed usage run `ml_predict.R -h`

    usage: /home/springer/zhoux379/git/nf/bin/mmm/ml_predict.R [-h] fm fi fo
    
    Make predictions using trained model on given dataset
    
    positional arguments:
      fm          (ML) model file
      fi          input dataset
      fo          output file to save predictions
    
    optional arguments:
      -h, --help  show this help message and exit

