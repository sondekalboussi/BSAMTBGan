# BSAMTBGan (Belgian SouthAfrican MycoBacterium Tuberculosis Genome annotation)
Python pipeline for MTB reference based genome assembly

# Pipeline workflow
    Trimming the bad quality reads using Trimmomatic

    Mapping using bwa and novocraft 

    Realignment around the ondels and base quality qualibration using GATK

    Mark and remove PCR duplicates using picard

    Merge bam files in case we have one sample having reads coming from different libraries

    Check genome coverage and reads mappability: we accept reads with genome coverage >=40 and/or reads mappability >=90

    Joint varant call using GATK haplotypeCaller

    Structural variant call using Delly

    Gene annotation using Annovar

Converting the gtf file into genepred file for Annovar gene annotation must be run on linux OS (no mac/windows version),please install libpng12-0 :sudo apt-get install to be able to run gtfToGenePred code (http://hgdown
load.soe.ucsc.edu/admin/exe/linux.x86_64/).


You need to provide a table that contains columns in the following order :file path, ID, PL, SM and LB corresponding to identifiers needed to add readgroups to run trimmomatic (Please check Data folder to find a template of the table).









Please cite these scripts as follows: S. Kalboussi, A. Dippenaar & T.H. Heupink,(2017),BSAMTBGan.py https://github.com/sondekalboussi/BSAMTBGan
