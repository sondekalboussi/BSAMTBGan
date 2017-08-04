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









Please cite these scripts as follows: S. Kalboussi, A. Dippenaar & T.H. Heupink,(2017),BSAMTBGan.py https://github.com/sondekalboussi/BSAMTBGan
