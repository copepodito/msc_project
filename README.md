# MSc Bioinformatics with Systems Biology

This repository contains all code I used to manipulate and analyse the datasets of "The Great Antibody repertoire" by Briney *et. al* for my Master's Project (MP). Their work is publicly available [here](https://github.com/briney/grp_paper). Only annotated datasets for eight of the ten subjects originally considered in the study were used for my work, namely those present in the following list:
- 316188
  - Annotated data: [consensus CSVs](http://burtonlab.s3.amazonaws.com/sequencing-data/hiseq_2016-supplement/316188_HNCHNBCXY_consensus_UID18-cdr3nt-90_minimal_071817.tar.gz)
- 326650
  - Annotated data: [consensus CSVs](http://burtonlab.s3.amazonaws.com/sequencing-data/hiseq_2016-supplement/326650_HCGCYBCXY_consensus_UID18-cdr3nt-90_minimal_071817.tar.gz)
- 326737
  - Annotated data: [consensus CSVs](http://burtonlab.s3.amazonaws.com/sequencing-data/hiseq_2016-supplement/326737_HNKVKBCXY_consensus_UID18-cdr3nt-90_minimal_071817.tar.gz)
- 326780
  - Annotated data: [consensus CSVs](http://burtonlab.s3.amazonaws.com/sequencing-data/hiseq_2016-supplement/326780_HLH7KBCXY_consensus_UID18-cdr3nt-90_minimal_071817.tar.gz)
- 326797
  - Annotated data: [consensus CSVs](http://burtonlab.s3.amazonaws.com/sequencing-data/hiseq_2016-supplement/326797_HCGNLBCXY%2BHJLN5BCXY_consensus_UID18-cdr3nt-90_minimal_071817.tar.gz)
- 326907
  - Annotated data: [consensus CSVs](http://burtonlab.s3.amazonaws.com/sequencing-data/hiseq_2016-supplement/326907_HLT33BCXY_consensus_UID18-cdr3nt-90_minimal_071817.tar.gz)
- 327059
  - Annotated data: [consensus CSVs](http://burtonlab.s3.amazonaws.com/sequencing-data/hiseq_2016-supplement/327059_HCGTCBCXY_consensus_UID18-cdr3nt-90_minimal_071817.tar.gz)
- D103
  - Annotated data: [consensus CSVs](http://burtonlab.s3.amazonaws.com/sequencing-data/hiseq_2016-supplement/D103_HCGCLBCXY_consensus_UID18-cdr3nt-90_minimal_071817.tar.gz)
  
Scripts in this repository are implemented for eight directories, each one named after the number of the subject (e.g. `316188`) and containing their respective consensus CSVs. These directories should be placed within a parent directory named `~/Desktop/Project` for the original code to work. These requirements can be easily modified if desired (see commented code).
