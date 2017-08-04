# BSAMTBGan (Belgian SouthAfrican MycoBacterium Tuberculosis Genome annotation)
Python pipeline for MTB reference based genome assembly

# Pipeline workflow

   MTB strain identification using Spotyping 
    
   Trimming of the bad quality reads using Trimmomatic

   Mapping reads using bwa and novocraft to the H37Rv reference genome

   Realignment around the indels and base quality qualibration using GATK

    ark and remove PCR duplicates using picard

   Merge bam files in case we have one sample with reads coming from different libraries

   Check genome coverage and reads mappability: Only accept reads with genome coverage >=40 and/or reads mappability >=90

   Joint varant call using GATK haplotypeCaller

   Structural variant call using Delly
    
   Filter out all SNPs that arose in repetitive regions of the genome (please check Exclude.txt in Data directory)

   Gene annotation using Annovar

Converting the gtf file into genepred file for Annovar gene annotation must be run on linux OS (no mac/windows version),please install libpng12-0 :sudo apt-get install to be able to run  code
You can download the gtfToGenePred code from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/.


Please provide samples fastq files csv file that contains columns in the following order :Fastq file path, ID, PL, SM and LB for each sample corresponding to identifiers needed to add readgroups to run trimmomatic (Please check Data folder to find a template of the table).









Please cite the code as follows: S. Kalboussi, A. Dippenaar & T.H. Heupink,(2017),BSAMTBGan https://github.com/sondekalboussi/BSAMTBGan
