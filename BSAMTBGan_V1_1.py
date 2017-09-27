
"""Please make sure to install java and perl 
converting the gtf file into genepred file must be run on linux OS (no mac/windows version)
self.fastQFileData is a table that contains columns in the following order :fastq file path, ID, PL, SM and LB corresponding identifiers needed to add readgroups to bam files
install libpng12-0 :sudo apt-get install
Change the genome fasta header to NC-000962-3-H37Rv, the olde header will cause a bug for novoalign mapper
change contigs name in the dbSNP to  NC-000962-3-H37Rv
remove trailing white space at dbsnp.csv fiel
"""
import os
import sys
import csv
from itertools import product
import fileinput
from datetime import datetime
start_time = datetime.now()
new=open("/Users/sondeskalboussi/Desktop/Time1.txt","w")
class BSMTBGan(object):
    def __init__(self,global_Dir,user_ref,user_ref_fasta,File_table,mappers):
      
      #Settings parameters and paths
        self.global_Dir=global_Dir# path to home directory in the cluster or local machine
        self.results=self.global_Dir+"/Results/"#results folder
        self.tools="/Users/sondeskalboussi/Tools"#Tools binaries in the cluster
        self.ref_gen_Dir="/Users/sondeskalboussi/Desktop/BSATBGan_pipeline/Reference/"#H37rv reference
        self.genome=user_ref
        self.fasta=user_ref_fasta
        self.ref_genome="/Users/sondeskalboussi/Desktop/BSATBGan_pipeline/Reference/"+self.genome+"/FASTA/"#fasta file of reference genome
        self.ref="/Users/sondeskalboussi/Desktop/BSATBGan_pipeline/Reference/"+self.genome+"/FASTA/"+self.fasta#fasta file of reference genome
        self.novocraft=self.tools+"/novocraft/"# path to the novocraft software
        self.bwa=self.tools+"/bwa-0.7.15/"#path to bwa mem software
        self.trimmomatic=self.tools+"/Trimmomatic-0.36/"
        self.GATK=self.tools+"/GenomeAnalysisTK.jar"
        self.samtools="/samtools/"
        self.picard=self.tools+"/picard.jar"
        self.bedtools="bedtools/"
        self.bcftools=self.tools+"/bcftools-1.4/"
        self.Delly=self.tools+"/delly/src/delly"
        self.annovar=self.tools+"/annovar/"
        self.mappers=mappers
        self.dbSNP="/Users/sondeskalboussi/Desktop/BSATBGan_pipeline/Reference/dbSNP/dbSNP.vcf"
        self.illumina_adapters="/Users/sondeskalboussi/Desktop/BSATBGan_pipeline/illumina_adapters.fna.fasta"
        self.fastQFileData=File_table
        self.gtfToGenePred="/Users/sondeskalboussi/Desktop/BSATBGan_pipeline/gtfToGenePred"#bin file
        self.H37RV_gtf="/Users/sondeskalboussi/Desktop/BSATBGan_pipeline/Mycobacterium_tuberculosis_h37rv.ASM19595v2.36.gtf"
        self.TBdb="/Users/sondeskalboussi/Desktop/BSATBGan_pipeline/dbMTB"
        self.spotyping=self.tools+"/SpoTyping-2.1/SpoTyping-v2.1-commandLine/"
        self.Sample=[]
        self.libraryPE={}
        self.librarySE={}
        self.valide_bam=[]#bam files validated by statistics for variant call
    #self.spo_outp=spo_outp#spoligotyping output directory is the working directory by default, if u put the full path u ll have an error message so you need to


    print 
    print
    print "_____________________________________________________________________________________________"
    print "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
    print "BSMTBGan - Belgian SouthAfrican Tuberculosis genome Annotation pipeline                   "
    print "V1.0, Created by Sondes kalboussi, May 2017, email: sondesbioinf@gmail.com         "
    print                       "Copyright Sondes kalboussi 2017 "
    print "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
    print "_____________________________________________________________________________________________"
    print ""
    
    #Process the reference genome for the later use for mapping andGATK_
    def process_ref_genome(self):
           print ("""
                ============================================================================
                     Reference genomes processing for the annotation is starting"!
                ============================================================================
            """)
           for fl in os.listdir(self.ref_genome):
               if fl.endswith(".fa") or fl.endswith(".fasta"):
                        os.chdir(self.ref_genome)
                        ID=os.path.splitext(fl)[0]
                        os.system("{}bwa index {}""".format(self.bwa,fl))#index the reference fasta for bwa
                        os.system("{}novoindex {}.nix {} ".format(self.novocraft,ID,fl))#index the reference fasta for novoalign
                        os.system("samtools faidx {}""".format(fl))
                        os.system("java -jar {} CreateSequenceDictionary R={} O={}.dict""".format(self.picard,fl,ID))#GATK
           return "Processed genome files are ready !"
    # Strain identification using Spotyping
    def spoligotyping(self):
        #use Spotyping software to classify Mtb strains  
        file=csv.reader(open(self.fastQFileData,'rb'))
        next(file, None)
        print ("""
                  ========================================================
                             Strain identification starts now!
                  ========================================================
        """)
        if not os.path.isdir(self.results+"Spoligotyping"):
            os.makedirs(self.results+"Spoligotyping")
        for row in file:
                SM=row[3]
                LB=row[4]
                if "R1" not in row[0] and "R2" not in row[0]:
                        readSE=row[0]
                        output_spolig_SE="Desktop/BSAMGan_test1/Results/Spoligotyping/"+SM+"_"+LB+"_SE_spo.txt"
                        os.system("""python {}SpoTyping.py {} --output {}""".format(self.spotyping,readSE,output_spolig_SE))
                if "R1" or "R2" in row[0] and SM in row[0] :
                        readPE=row[0]
                        if "R1" in readPE:
                            read1=readPE
                        if "R2"in readPE:
                            read2=readPE
                            output_spolig_PE="Desktop/BSAMGan_test1/Results/Spoligotyping/"+SM+"_"+LB+"_PE_spo.txt"
                            os.system("""python {}SpoTyping.py {} {} --output {}""".format(self.spotyping,read1,read2,output_spolig_PE))
        return "Spoligotyping is finished!"
    #Reads Trimming 
    def trimming(self):
        file=csv.reader(open(self.fastQFileData,'rb'))
        next(file, None)
        print ("""
                  ====================================
                       The trimming starts now!
                  ====================================
        """)

        if not os.path.isdir(self.results+"Trimming"):
               os.makedirs(self.results+"Trimming")
        for row in file:
                SM=row[3]
                if SM not in self.Sample:
                    self.Sample.append(SM)
                
                ID=row[1]
                PL=row[2]
                LB=row[4]
                if "R1" not in row[0] and "R2" not in row[0]:
                    readSE=row[0]
                    output_trim=self.results+"Trimming/"+SM+"_"+ID+"_"+LB+"_"+PL+"_SE_trimed.fq"
                    os.system("""java -jar {}trimmomatic-0.36.jar SE -threads 4 -phred33 {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 MINLEN:36""".format(self.trimmomatic,readSE,output_trim))
                if "R1" or "R2" in row[0] and SM in row[0] :
                    read=row[0]
                    if "R1" in read:
                        read1=read
                        LB1=row[4]
                    if "R2"in read:
                        read2=read
                        LB2=row[4]
                        if LB1==LB2:
                            output1=self.results+"Trimming/"+SM+"_"+ID+"_"+LB+"_"+PL+"_R1.PE_paired_trimed.fq"
                            output2=self.results+"Trimming/"+SM+"_"+ID+"_"+LB+"_"+PL+"_R1.PE_unpaired_trimed.fq"
                            output3=self.results+"Trimming/"+SM+"_"+ID+"_"+LB+"_"+PL+"_R2.PE_paired_trimed.fq"
                            output4=self.results+"Trimming/"+SM+"_"+ID+"_"+LB+"_"+PL+"_R2.PE_unpaired_trimed.fq"
                            os.system("""java -jar {}trimmomatic-0.36.jar PE -threads 4 -phred33 {} {} {} {} {} {} ILLUMINACLIP:TruSeq3-PE:2:30:10 LEADING:3 TRAILING:3 MINLEN:36""".format(self.trimmomatic,read1,read2,output1,output2,output3,output4))#change extension of trimmomatic for the cluster
    # Create a dictionary where the keys are samples ID and value is the reads ID when they are coming from more than sequencing library so later you can merge all the output for the same sample 
        for fl in os.listdir(self.results+"Trimming/"):
                os.chdir(self.results+"Trimming/")
                if fl.endswith("R1.PE_paired_trimed.fq"):
                                    IDPE=fl.split("R1.PE_paired_trimed.fq")[0]
                                    sm=IDPE.split("_")[0]
                                    lib=IDPE.split("_")[2]
                                    if sm not in self.libraryPE.keys():
                                        self.libraryPE[sm]=[lib]
                                    else:
                                        self.libraryPE[sm].append(lib)
                if fl.endswith("_SE_trimed.fq"):
                                    IDSE=fl.split("_SE_trimed.fq")[0]
                                    sm=IDSE.split("_")[0]
                                    lib=IDSE.split("_")[2]
                                    if sm not in self.librarySE.keys():
                                        self.librarySE[sm]=[lib]
                                    else:
                                        self.librarySE[sm].append(lib)
        
        return "Trimming is done!"
   
   # Mapping reads using BWA-mem and Novoalign, mappers is a list of software that we want to use 
    def Mapping(self):
        print ("""
              ========================
                Mapping starts now!
              ========================
            """)
        for fl in os.listdir(self.results+"Trimming/"):
          os.chdir(self.results+"Trimming/")
          if ".fq" in fl:
            SM=fl.split("_")[0]
            ID=fl.split("_")[1]
            PL=fl.split("_")[3]
            LB=fl.split("_")[2]
            readGroup="@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL
            input=os.path.join(self.results+"Trimming/",fl)
            for map in self.mappers:# mappers is a list of mapping software name (BWA-mem and novoalign in our study)
                    if not os.path.isdir(self.results+map+"/Alignment/Final_Bam"):
                        os.makedirs(self.results+map+"/Alignment/Final_Bam")
                    #Mapping the SE reads
                    if "SE_trimed" in input:
                        output_map=self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".sam"
                        if map=="BWA":
                            os.system("""{}bwa mem -M -t 4 -R "{}" {} {} > {}""".format(self.bwa,readGroup,self.ref,input,output_map))
                            os.system("""samtools view -Sb -o {} {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_PE_"+map+".bam",output_map))
                        if map=="Novoalign":
                            os.system("""{}novoalign -d {} -f {} -o SAM '{}' > {}""".format(self.novocraft,self.ref_genome+self.fasta.split(".")[0]+".nix",input,readGroup,output_map))
                            os.system("""samtools view -Sb -o {} {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam",output_map))
                    #Mapping the PE reads
                    if "PE_paired_trimed" in input:
                        if "R1" in input:
                            R1=input
                        if "R2" in input:
                            R2=input
                            output_map=self.results+map+"/Alignment/"+SM+"_"+LB+"_PE_"+map+".sam"
                            if map=="BWA":
                                os.system("""{}bwa mem -M -t 4 -R "{}" {} {} {} > {}""".format(self.bwa,readGroup,self.ref,R1,R2,output_map))
                                os.system("""samtools view -Sb -o {} {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_PE_"+map+".bam",output_map))
                            if map=="Novoalign":
                              os.system("""{}novoalign -d {} -f {} {} -o SAM '{}' > {}""".format(self.novocraft,self.ref_genome+self.fasta.split(".")[0]+".nix",R1,R2,readGroup,output_map))
                              os.system("""samtools view -Sb -o {} {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_PE_"+map+".bam",output_map))
        for map in self.mappers:
                for fl in os.listdir(self.results+map+"/Alignment"):
                    os.chdir(self.results+map+"/Alignment")
                    if fl.endswith(".sam"):
                            os.system("""JAVA -jar {} SortSam INPUT={} OUTPUT={} SORT_ORDER=coordinate""".format(self.picard,fl,fl.replace(".sam","_sorted.bam")))
                            os.system("""JAVA -jar {} BuildBamIndex INPUT={}""".format(self.picard,fl.replace(".sam","_sorted.bam")))
        return "Mapping is done"
    #Mark and remove the PCR duplicates
    def PCR_dup_mark(self):
           
            print ("""
                            ==========================================
                                PCR duplication marking begins!
                            ==========================================
                          """)
            for map in self.mappers: 
                for fl in os.listdir(self.results+map+"/Alignment/"):
                    os.chdir(self.results+map+"/Alignment/")
                    if fl.endswith(map+"_sorted.bam"):
                            input=os.path.join(self.results+map+"/Alignment/",fl)
                            output=os.path.join(self.results+map+"/Alignment/",fl.replace(map+"_sorted.bam",map+"_dedup.bam"))
                            metric=os.path.join(self.results+map+"/Alignment/",fl.replace(map+"_sorted.bam",map+"_dedup.bam.metrics"))
                            os.system("""java -jar {} MarkDuplicates INPUT={} OUTPUT={} METRICS_FILE={}""".format(self.picard,input,output,metric))#add{}
                            os.system("""java -jar {} BuildBamIndex INPUT={}""".format(self.picard,output))
            return "PCR duplicates removal is done!"
                        
   # Realign around the indels using GATK the final output is a sorted indexed bam file
    def realignment(self):                 
        print ("""
                =====================================================
                      Realignment around the indels using GATK
                =====================================================
            """)
        for map in self.mappers:           
            for fl in os.listdir(self.results+map+"/Alignment/"):
                  os.chdir(self.results+map+"/Alignment/")
                  if fl.endswith("_dedup.bam"):
                        input=os.path.join(self.results+map+"/Alignment/",fl)
                        output1=os.path.join(self.results+map+"/Alignment/",fl.replace("_dedup.bam",".intervals"))
                        output2=os.path.join(self.results+map+"/Alignment/",fl.replace("_dedup.bam","_realg.bam"))
                        output3=os.path.join(self.results+map+"/Alignment/",fl.replace("_dedup.bam","_realg_sorted.bam"))
                        os.system("""java -jar {} -T RealignerTargetCreator -nt 4 -R {} -I {} -o {}""".format(self.GATK,self.ref,input,output1))#change thread
                        os.system("""java -jar {} -T IndelRealigner -R {} -I {} -targetIntervals {} -o {}""".format(self.GATK,self.ref,input,output1,output2))
                        os.system("""JAVA -jar {} SortSam INPUT={} OUTPUT={} SORT_ORDER=coordinate""".format(self.picard,output2,output3))
                        #os.system("""samtools sort -o {} {}""".format(output3,output2))
                        os.system("""java -jar {} BuildBamIndex INPUT={}""".format(self.picard,output3))
        return "realignment around indels is done!"
    
    #Base quality score recalibration
    def base_qual_recal(self):                   
        for map in self.mappers:
             print ("""
                =================================================================
                    {} alignment base quality score recalibration starts now!
                =================================================================
             """.format(map))                     
             for fl in os.listdir(self.results+map+"/Alignment/"):
                  os.chdir(self.results+map+"/Alignment/")
                  if fl.endswith(map+"_realg_sorted.bam"):
                        input=os.path.join(self.results+map+"/Alignment/",fl)
                        output1=os.path.join(self.results+map+"/Alignment/",fl.replace("_realg_sorted.bam","_recal_data.table"))
                        output2=os.path.join(self.results+map+"/Alignment/",fl.replace("_realg_sorted.bam","_recal.bam"))
                        output3=os.path.join(self.results+map+"/Alignment/",fl.replace("_realg_sorted.bam","_recal_sorted.bam"))
                        os.system("""java -jar {} -T BaseRecalibrator -nct 4 -R {} -I {} -knownSites {} -o {}""".format(self.GATK,self.ref,input,self.dbSNP,output1))#add path to GATK cluster
                        os.system("""java -jar {} -T PrintReads -nct 4 -R {} -I {} -BQSR {}  -o {}""".format(self.GATK,self.ref,input,output1,output2))
                        #os.system("""samtools sort -o {} {}""".format(output3,output2))
                        os.system("""JAVA -jar {} SortSam INPUT={} OUTPUT={} SORT_ORDER=coordinate""".format(self.picard,output2,output3))
                        os.system("""java -jar {} BuildBamIndex INPUT={}""".format(self.picard,output3))

        return "Base recalibration is finished!"
       #Check if there is samples from more than one sequencing library if yes merge them and move them to Final_bam folder otherwise use bam file in Alignment directory with the extension final_mapper.bam
    def merge_bam(self):
        for map in self.mappers:                
            #samples that do need merging are moved to a new final folder if sample id has one lib value
                for key,value in self.libraryPE.iteritems():
                    if len(value)==1:#no merging just move bam file to final_Bam folder
                        for i in value:
                            os.system("cp "+self.results+map+"/Alignment/"+key+"_"+i+"_PE_"+map+"_recal_sorted.bam "+self.results+map+"/Alignment/Final_Bam/"+key+"_"+i+"_PE_"+map+"_recal_sorted.bam")
                    else: #merge all the bams into final one
                         print ("""            
                      ====================================
                       Merging PE Bam files for {} mapper
                      ====================================
                             """.format(map))
                         multi_libraryPE=[ ]
                         to_be_mergedPE=""
                         for i in value:
                            multi_libraryPE.append(self.results+map+"/Alignment/"+key+"_"+i+"_PE_"+map+"_recal_sorted.bam")
                                        #when i find the id in the file name write the file name multi
                         for i in multi_libraryPE:
                            to_be_mergedPE+=i+" "
                         output=self.results+map+"/Alignment/Final_Bam/"+key+"_PE_merged_"+map+".bam"
                         cmdPE="""samtools merge {} {}""".format(output,to_be_mergedPE)#add {}samtools
                         os.system(cmdPE)
                for key,value in self.librarySE.iteritems():
                    if len(value)==1:
                        for i in value:
                            os.system("cp "+self.results+map+"/Alignment/"+key+"_"+i+"_SE_"+map+"_recal_sorted.bam "+self.results+map+"/Alignment/Final_Bam/"+key+"_"+i+"_SE_"+map+"_recal_sorted.bam")
                    else:
                        print ("""            
                      =============================================                       
                           Merging SE Bam files for {} mapper
                      =============================================
                            """.format(map))
                        multi_librarySE=[ ]
                        to_be_mergedSE=""
                        for i in value:
                            multi_librarySE.append(self.results+map+"/Alignment/"+key+"_"+i+"_SE_"+map+"_recal_sorted.bam")
                                        #when i find the id in the file name write the file name multi
                        for i in multi_librarySE:
                            to_be_mergedSE+=i+" "
                        output=self.results+map+"/Alignment/Final_Bam/"+key+"_SE_merged_"+map+".bam"
                        cmdSE="""samtools merge {} {}""".format(output,to_be_mergedSE)
                        os.sys(cmdSE)
                 #sort and index the final bam                        
                for fl in os.listdir(self.results+map+"/Alignment/Final_Bam/"):
                            os.chdir(self.results+map+"/Alignment/Final_Bam/")
                            os.system("""java -jar {} BuildBamIndex INPUT={}""".format(self.picard,fl))
        return "The Fianl bam files are ready for further processing!"
      # Check the mapping quality only samples with >=90% mapped reads and/or genome coverage>=40% are accepted
    def mapping_stat(self):
      
            print ("""
                 ==================================================================
                   Genome coverage and reads mappability statistics are starting !
                 ==================================================================
                """)
            for map in self.mappers:
                MAP={}
                COV={}
                if not os.path.isdir(self.results+map+"/statistics"):
                    os.makedirs(self.results+map+"/statistics")
                if len(os.listdir(self.results+map+"/Alignment/Final_Bam/"))== 0:
                    Final_Bam=self.results+map+"/Alignment/"
                else:
                    Final_Bam=self.results+map+"/Alignment/Final_Bam/"
                for fl in os.listdir(Final_Bam):
                        os.chdir(Final_Bam)
                        if fl.endswith("_recal_sorted.bam") or "_merged_" in fl:
                            input=os.path.join(Final_Bam,fl)
                            output1=self.results+map+"/statistics/"+fl.replace(".bam",".coverage")
                            output2=self.results+map+"/statistics/"+fl.replace(".bam",".flagstat")
                            output3=self.results+map+"/statistics/"+fl.replace(".bam",".stats")
                            os.system("""java -jar {} -T DepthOfCoverage -R {} -I {} -o {} --omitDepthOutputAtEachBase --omitIntervalStatistics --omitLocusTable""".format(self.GATK,self.ref,input,output1))
                            os.system("""samtools flagstat {} > {}""".format(input,output2))
                            os.system("""samtools stats {} > {}""".format(input,output3))
                for fl in os.listdir(self.results+map+"/statistics/"):
                   os.chdir(self.results+map+"/statistics/")
                   ID=fl.split("_")[0]
                   if fl.endswith(".flagstat"):
                       for line in open(fl).readlines():
                            line=line.strip()
                            if "mapped (" in line:
                                x=line.split("mapped")[1].split(":")[0].split("(")[1].split("%")[0]
                                y=float(x)
                                MAP[ID]=y
                   if fl.endswith("_summary"):
                       for line in open(fl).readlines()[1:2]:
                                line=line.strip()
                                x=line.split("\t")[2]
                                y=float(x)
                                COV[ID]=y
                new=open(self.results+map+"/statistics/mapping_stat.txt","w")
                new1=open(self.results+map+"/statistics/Failed_mapping_stat.txt","w")
                new.write("""{}{}{}{}{}{}""".format("Sample","\t","%Mapped_reads","\t","Genome_coverage mean", "\n"))
                new1.write("""{}{}{}{}{}{}""".format("Sample","\t","%Mapped_reads","\t","Genome_coverage mean", "\n"))
                for key in MAP.keys():
                    if key in COV.keys():
                        if MAP[key]>=90 and COV[key]>=40 :
                            if key not in self.valide_bam:
                                self.valide_bam.append(key)
                            new=open(self.results+map+"/statistics/mapping_stat.txt","a")
                            new.write("""{} {}  {}{}""".format(key,MAP[key],COV[key],"\n"))
                    
                        else:
                            new1=open(self.results+map+"/statistics/Failed_mapping_stat.txt","a")
                            new1.write("""{} {}  {}{}""".format(key,MAP[key],COV[key],"\n"))
                new.close()
                new1.close()
            return " Statistics are done!"
    
    # Joint variant call SNP and indels using GATK
    def joint_variant_calling(self):
            print ("""
                  ==================================================================
                    Joint variant calling (SNP and indels) using GATK is starting !
                  ==================================================================
            """)                        
            for map in self.mappers:
                if not os.path.isdir(self.results+map+"/Joint_Variants"):
                    os.makedirs(self.results+map+"/Joint_Variants")
                if len(os.listdir(self.results+map+"/Alignment/Final_Bam/"))== 0:
                    Final_Bam=self.results+map+"/Alignment/"
                else:
                    Final_Bam=self.results+map+"/Alignment/Final_Bam/"
                for i in self.valide_bam:
                    for fl in os.listdir(Final_Bam):
                        if i in fl:
                            if fl.endswith("_recal_sorted.bam") or "_merged_" in fl:
                                os.chdir(Final_Bam)
                                input=os.path.join(Final_Bam,fl)
                                output=os.path.join(self.results+map+"/Joint_variants/",fl.replace(".bam","_GATK_snps_indels.vcf"))
                                os.system("""java -jar {} -T HaplotypeCaller -R {} -I {} -stand_call_conf 30 -o {}""".format(self.GATK,self.ref,input,output))#add path gatk cluster
            return "Joint variant calling (SNP and indels) using GATK is done!"
    #non-model organism, so variant re calibration isn't possible
    def joint_variant_calling_hard_filtering(self):
         print ("""
            ==================================================================
                  Joint_variant_calling_hard_filtering is starting !
            ==================================================================
            """)
         for map in self.mappers:
             for fl in os.listdir(self.results+map+"/Joint_Variants"):
                 os.chdir(self.results+map+"/Joint_Variants")
                 id=fl.split("_")[0]
                 if fl.endswith(".vcf"):
                     input=os.path.join(self.results+map+"/Joint_Variants",fl)
                     output1=self.results+map+"/Joint_variants/"+id+"_raw_snps.vcf"
                     output11=self.results+map+"/Joint_variants/"+id+"_filtered_snps.vcf"
                     output2=self.results+map+"/Joint_variants/"+id+"_raw_indels.vcf"
                     output22=self.results+map+"/Joint_variants/"+id+"_filtered_indels.vcf"
                     #variant selection
                     os.system("""java -jar {} -T SelectVariants -R {} -V {} -selectType SNP -o {}""".format(self.GATK,self.ref,input,output1))
                     os.system("""java -jar {} -T SelectVariants -R {} -V {} -selectType INDEL -o {}""".format(self.GATK,self.ref,input,output2))
                     #variant filtration
                     os.system("""java -jar {} -T VariantFiltration -R {} -V {} --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o {}""".format(self.GATK,self.ref,output1,output11))
                     os.system("""java -jar {} -T VariantFiltration -R {} -V {} --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter" -o {}""".format(self.GATK,self.ref,output2,output22))
         return " variant_calling_hard_filtering is done!"
    #number of SNPs identified with average quality scores and average mapping quality for filtered snp.csv and raw.vcf, parse the MQ mapping quality, QUAl field: quality score, wc - l nbre of snp then compute average
    def SNP_statistics(self):
        print ("""
            ==================================================================
                  SNP statistics is starting !
            ==================================================================
            """)
        for map in self.mappers:
            new=open(self.results+map+"/statistics/VCF_stat.txt","w")
            new.write("""{}{}{}{}{}{}{}{}""".format("Sample","\t","Nbre_SNP","\t","AV_MQ","\t","AV_QUAL","\n"))
            for fl in os.listdir(self.results+map+"/Joint_variants/"):
                os.chdir(self.results+map+"/Joint_variants/")
                f=os.path.join(self.results+map+"/Joint_variants/",fl)
                if f.endswith("_filtered_snps.vcf"):
                    Nbre_SNP=0
                    AV_MQ=0
                    AV_QUAL=0
                    for lines in filter(None, open(f).readlines()):
                        if not lines.startswith("##") and not lines.startswith("#") and "my_snp_filter" not in lines:
                            Nbre_SNP+=1
                            line=lines.split("\t")
                            AV_QUAL+=float(line[5])/Nbre_SNP
                            L=line[7].split(";")[8]
                            if "MQ" in L:
                                AV_MQ+=float(line[7].split(";")[8].split("=")[1])/Nbre_SNP
                    new.write("""{}{}{}{}{}{}{}{}""".format(fl.split("_")[0],"\t",Nbre_SNP,"\t",AV_MQ,"\t",AV_QUAL,"\n"))
            new.close()
        return "VCF statistics is done"
    #structural variant identification(big deletion) using DELLY
    def Genotype_structural_variation_calling(self):
        print ("""
                  ======================================================
                    Genotupe structural_variation_calling is starting !
                  ======================================================
        """)                    
        for map in self.mappers:
            if not os.path.exists(self.results+map+"/Struc_Variants"):
                    os.makedirs(self.results+map+"/Struc_Variants")
            if len(os.listdir(self.results+map+"/Alignment/Final_Bam/"))== 0:
                Final_Bam=self.results+map+"/Alignment/"
            else:
                Final_Bam=self.results+map+"/Alignment/Final_Bam/"
            for i in self.valide_bam:
                for fl in os.listdir(Final_Bam):
                    if i in fl:
                        if fl.endswith("_recal_sorted.bam") or "_merged_" in fl:
                            os.chdir(Final_Bam)
                            input=input=os.path.join(Final_Bam,fl)
                            output=self.results+map+"/Struc_Variants/"+fl.replace("_realg_sorted.bam","_SV")
                            os.system("""delly call -t DEL -g {} -o {}.bcf {}""".format(self.ref,output,input))#add delly
                            os.system("""bcftools view {}.bcf > {}.vcf """.format(output,output))#convert bcf to vcf
        return "GSV is done!"        
    #SNPs and Indels are annotated as intergenic or genic and for amino acid level changes using ANNOVAR
    def Gene_annotation_annovar(self):
        """Create the gene annotation databases:gtf file and genome assembly fasta
           convert the VCF into annovar format
           """
        print ("""
              ====================================
                   Annotation is starting !
              ====================================
        """)
        '''
        if not os.path.isdir(self.TBdb):
            os.makedirs(self.TBdb)
        ##download the genome assembly file from Ensembl
        os.system("wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.toplevel.fa.gz -P "+self.TBdb)
        os.system("gunzip "+self.TBdb+"/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.toplevel.fa.gz")
        #convert gene definition file into gene prediction file
        os.system("""chmod +x {}""".format(self.gtfToGenePred))
        for fl in os.listdir(self.TBdb):
            if not os.path.exists(self.TBdb+"/H37Rv_refGene.txt")  and not os.path.exists(self.TBdb+"/H37Rv_refGeneMrna.fa"):
                os.system(""".{} -genePredExt {} {}""".format(self.gtfToGenePred,self.H37RV_gtf,self.TBdb+"/H37Rv_refGene.txt")) #run a bin file (first check if you have lib2 installed, only on linux OS
        ##Generate a transcript fasta file from the genome assembly file
                os.system("""perl {}/retrieve_seq_from_fasta.pl --format refGene --seqfile {} {} --out {}""".format(self.annovar,self.ref,self.TBdb+"/H37Rv_refGene.txt",self.TBdb+"/H37Rv_refGeneMrna.fa"))'''
        ##run annovar gene annotation()/ Convert the VCF to annovar input format file if you want to use annovar command with avinput but u need to modify the command first:
        for map in self.mappers:
            if not os.path.isdir(self.results+map+"/Annotation"):
                os.makedirs(self.results+map+"/Annotation")
            for fl in os.listdir(self.results+map+"/Joint_Variants/"):
                os.chdir(self.results+map+"/Joint_Variants/")
                if fl.endswith("_snps_indels.vcf"):
                        input=fl
                        id=fl.split("_")[0]
                        #output=self.results+map+"/Joint_Variants/"+fl.replace("_snps_indels.vcf",".avinput")
                        #os.system("""perl {}convert2annovar.pl -format vcf4 {} > {}""".format(self.annovar,input,output))
                        #os.system("""sed -i -e 's/NC-000962-3-H37Rv/Chromosome/g' {}""" .format(input))'''
                        #os.system("""perl {}table_annovar.pl {} {} --vcfinput --outfile {} -buildver H37Rv --protocol refGene -operation g """.format(self.annovar,input,self.TBdb,id+"_"+map))
                        #os.system("""mv {} {}""".format(self.results+map+"/Joint_Variants/"+id+"_"+map+".H37Rv_multianno.vcf",self.results+map+"/Annotation/"+id+"_"+map+".H37Rv_multianno.vcf"))
                        os.system("""perl {}table_annovar.pl {} {} --vcfinput --outfile {} -buildver H37Rv --protocol refGene -operation g """.format(self.annovar,input,self.TBdb,id+"_"+map))
                        os.system("""mv {} {}""".format(self.results+map+"/Joint_Variants/"+id+"_"+map+".H37Rv_multianno.vcf",self.results+map+"/Annotation/"+id+"_"+map+".H37Rv_multianno.vcf"))
        return " Annovar Gene annotation is done!"
        
def main():
        pipeline=BSMTBGan("/Users/sondeskalboussi/Desktop/BSAMGan_test1","Mycobacterium_tuberculosis_h37rv","Mycobacterium_tuberculosis_h37rv.fa","/Users/sondeskalboussi/Desktop/BSAMGan_test/Fastq.csv", mappers=["BWA"])
        try:
            #print pipeline.process_ref_genome()
            #print pipeline.spoligotyping()
            #print pipeline.trimming()
            #print pipeline.Mapping()
            #print pipeline.PCR_dup_mark()
            #print pipeline.realignment()
            #print pipeline.base_qual_recal()
            #print pipeline.merge_bam()
            #print pipeline.mapping_stat()
            #print pipeline.joint_variant_calling()
            #print pipeline.joint_variant_calling_hard_filtering()
            #print pipeline.SNP_statistics()
            #print pipeline.Genotype_structural_variation_calling()
            #print pipeline.Gene_annotation_annovar()
        except IOError as e:
            print("I/O error: {0}".format(e))    
    
if __name__ == '__main__':main()
end_time = datetime.now()
new.write('BWA pipeline execution time: {}'.format(end_time - start_time))
new.close()

