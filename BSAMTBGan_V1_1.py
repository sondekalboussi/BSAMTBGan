
"""Please make sure to install java and perl 
converting the gtf file into genepred file must be run on linux OS (no mac/windows version)
self.fastQFileData is a table that contains columns in the following order :fastq file path, ID, PL, SM and LB corresponding identifiers needed to add readgroups to bam files
install libpng12-0 :sudo apt-get install 
"""
import os
import sys
import csv
from itertools import product
class BSMTBGan:
    def __init__(self,working_Dir,user_ref,user_ref_fasta,File_table,mappers):
      
      #Settings parameters and paths
          self.global_Dir=working_Dir# path to your working directory in the cluster
          self.results=self.global_Dir+"/Results/"#results folder
          self.tools="/opt/software/"#Tools binaries in the cluster
          self.ref_gen_Dir="/home/shared_data_ghi/CWGSMTB/Software/USAP/Reference/"#H37rv reference
          self.genome=user_ref
          self.fasta=user_ref_fasta
          self.ref_genome=self.ref_gen_Dir+self.genome+"/FASTA/"#fasta file of reference genome
          self.ref=self.ref_gen_Dir+self.genome+"/FASTA/"+self.fasta#fasta file of reference genome
          self.novocraft=self.tools+"/novocraft/3.02.07/bin/"# path to the novocraft software
          self.bwa=self.tools+"0.7.12/bin/"#path to bwa mem software
          self.trimmomatic=self.tools+"Trimmomatic-0.36/"
          self.GATK=self.tools+"gatk/GATK_3.6.0/"
          self.samtools=self.tools+"samtools/1.3.1/"
          self.picard=self.tools+"picard/picard-tools-2.3.0/"
          self.bedtools=self.tools+"BEDTools/2.16.2/"
          self.bcftools=self.tools+"bcftools/bcftools/1.4/"
          self.Delly=self.globalDir+"/Software/BSAMTBan/Tools_new/delly/src/"
          self.annovar=self.tools+"ANNOVAR/2017-06-01/"
          self.mappers=mappers
          self.dbSNp="/home/shared_data_ghi/CWGSMTB/Software/USAP/Reference/dbSNP/"
          self.illumina_adapters="/home/shared_data_ghi/CWGSMTB/Software/USAP/Tools/illumina_adapters.fna.fasta"
          self.fastQFileData=File_table
          self.gtfToGenePred="/gtfToGenePred/"#bin file
          self.H37RV_gtf=+"/Mycobacterium_tuberculosis_h37rv.ASM19595v2.36.gtf"
          self.TBdb=+"/dbMTB/"
          self.spotyping="/home/shared_data_ghi/CWGSMTB/Software/BSAMTBan/Tools/SpoTyping-2.1/SpoTyping-v2.1-commandLine/"
    print 
    print
    print "_____________________________________________________________________________________________"
    print "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
    print "BSMTBGan - Belgian SouthAfrican MTB genome Annotation pipeline                   "
    print "V1.0, Created by Sondes kalboussi, May 2017, email: sondes777@hotmail.com         "
    print "Copyright Sondes kalboussi 2017                                                           "
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
                        os.system("{}bwa index {}""".format(self.BWA,fl))#index the reference fasta for bwa
                        os.system("{}novoindex {}.nix {} ".format(self.novocraft,ID,fl))#index the reference fasta for novoalign
                        os.system("samtools faidx {}""".format(fl))
                        os.system("java -jar {} CreateSequenceDictionary R={}.fasta O={}.dict""".format(self.picard,fl,ID))#GATK
 
    # Strain identification using Spotyping
    def spoligotyping(self):
        global Sample
        Sample=[]
        #use Spotyping software to classify Mtb strains  
        file=csv.reader(open(self.global_Dir+"/"+self.fastQFileData,'rb'))
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
                if SM not in Sample:
                    Sample.append(SM)
                if "R1" not in row[0] and "R2" not in row[0]:
                        readSE=row[0]
                        output_spolig_SE="Spoligotyping/"+SM+"_"+LB+"_SE_spo.txt"
                        os.system("""python2.7 {}SpoTyping.py {} -O Results -o {}""".format(self.spotyping,readSE,output_spolig_SE))
                if "R1" or "R2" in row[0] and SM in row[0] :
                        readPE=row[0]
                        if "R1" in readPE:
                            read1=readPE
                        if "R2"in readPE:
                            read2=readPE
                            output_spolig_PE="Spoligotyping/"+SM+"_"+LB+"_PE_spo.txt"
                            os.system("""python2.7 {}SpoTyping.py {} {} -O Results -o {}""".format(self.spotyping,read1,read2,output_spolig_PE))
        return "Spoligotyping is finished!"
    #Reads Trimming 
    def trimming(self):
        global multi_libraryPE,multi_librarySE,Sample
        multi_libraryPE=dict()
        multi_librarySE=dict()
        Sample=list()
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
                if SM not in Sample:
                    Sample.append(SM)
                ID=row[1]
                PL=row[2]
                LB=row[4]
                if "R1" not in row[0] and "R2" not in row[0]:
                    readSE=row[0]
                    output_trim=self.results+"Trimming/"+SM+"_"+LB+"_SE_trimed.fq"
                    os.system("""java -jar {}trimmomatic-0.36.jar SE -threads 8 -phred33 {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(self.trimmomatic,readSE,output_trim))
                if "R1" or "R2" in row[0] and SM in row[0] :
                    read=row[0]
                    if "R1" in read:
                        read1=read
                        LB1=row[4]
                    if "R2"in read:
                        read2=read
                        LB2=row[4]
                        if LB1==LB2:
                            output1=self.results+"Trimming/"+SM+"_R1_"+LB1+"_PE_paired_trimed.fq"
                            output2=self.results+"Trimming/"+SM+"_R1_"+LB1+"_PE_unpaired_trimed.fq"
                            output3=self.results+"Trimming/"+SM+"_R2_"+LB1+"_PE_paired_trimed.fq"
                            output4=self.results+"Trimming/"+SM+"_R2_"+LB1+"_PE_unpaired_trimed.fq"
                            os.system("""java -jar {}trimmomatic-0.36.jar PE -threads 8 -phred33 {} {} {} {} {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(self.trimmomatic,read1,read2,output1,output2,output3,output4))#change extension of trimmomatic for the cluster
        
    # Create a dictionary where the keys are samples ID and value is the reads ID when they are coming from more than sequencing library so later you can merge all the output for the same sample 
        for fl in os.listdir(self.results+"Trimming/"):
                        os.chdir(self.results+"Trimming/")
                        if fl.endswith("_PE_paired_trimed.fq"): 
                                    IDPE=fl.split("_PE_paired_trimed.fq")[0]
                                    sm=IDPE.split("_")[0]
                                    lib=IDPE.split("_")[1]
                                    if sm not in multi_libraryPE.keys():
                                        multi_libraryPE[sm]=[lib]
                                    else:
                                        multi_libraryPE[sm].append(lib)
                        if fl.endswith("_SE_final_bwa.bam"):
                                    IDSE=fl.split("_SE_final_bwa.bam")[0]
                                    sm=IDSE.split("_")[0]
                                    lib=IDSE.split("_")[1]
                                    if sm not in multi_librarySE.keys():
                                        multi_librarySE[sm]=[lib]
                                    else:
                                        multi_librarySE[sm].append(lib)                  
        return "Trimming is done!"
   
   # Mapping reads using BWA-mem and Novoalign, mappers is a list of software that we want to use 
    def Mapping(self):
        print ("""
              ========================
                Mapping starts now!
              ========================
            """)
        file=csv.reader(open(self.fastQFileData,'rb'))
        next(file, None)
        for row in file:
            SM=row[3]
            ID=row[1]
            PL=row[2]
            LB=row[4]
            readGroup="@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL
        for map in self.mappers:# mappers is a list of mapping software name (BWA-mem and novoalign in our study)
            if not os.path.isdir(self.results+map+"/Alignment"): 
                os.makedirs(self.results+map+"/Alignment")
            for fl in os.listdir(self.results+"Trimming/"):
                  os.chdir(self.results+"Trimming/")
                  os.system("""gunzip -f -c {}""".format(fl))
                  if fl.endswith(".gz"):
                        input=os.path.join(self.results+"Trimming/",fl)
                    #Mapping the SE reads
                        if "SE" in input:
                            output_map=self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".sam"
                        if map=="BWA":
                            os.system("""{}bwa mem -M -t 8 -R "{}" {} {} > {}""".format(self.bwa,readGroup,self.ref,input,output_map))
                            os.system("""samtools view -Sb {} > {}""".format(output_map,self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam"))
                            os.system("""samtool sort {} > {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam",self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+"_sorted.bam"))#add{}/
                            os.system("""samtools index - > {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+"_sorted.bam",self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+"_sorted_idx.bam"))
                        if map=="Novoalign":
                            for fl in os.listdir(self.ref_genome):
                                            os.chdir(self.ref_genome)
                                            if fl.endswith(".nix"):  
                                                os.system("""{}novoalign -d {} -f {} -o SAM > {}""".format(self.novocraft,fl,input,output_map))
                                                os.system("""samtools view -Sb {} > {}""".format(output_map,self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam"))
                                                os.system("""samtool sort {} > {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam",self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+"_sorted.bam"))
                    #Mapping the PE reads
                        if "PE" in input:
                            if "R1" in input:
                                R1=input
                            if "R2" in input:
                                R2=input    
                                output_map=self.results+map+"/Alignment/"+SM+"_"+LB+"_PE_"+map+".sam"
                                if map=="BWA":
                                    os.system("""{}bwa mem -M -t 8 -R {} {} {} {} > {}""".format(self.bwa,readGroup,self.ref,R1,R2,output_map))
                                    os.system("""samtools view -Sb {} > {}""".format(output_map,self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam"))
                                    os.system("""samtool sort {} > {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam",self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+"_sorted.bam"))
                                if map=="Novoalign":
                                    for fl in os.listdir(self.ref_genome):
                                        os.chdir(self.ref_genome)
                                        if fl.endswith(".nix"):    
                                                os.system("""{}novoalign -d {} -f {} {} -i 200,50 -o SAM > {}""".format(self.novocraft,fl,R1,R2,output_map))
                                                os.system("""samtools view -Sb {} > {}""".format(output_map,self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam"))
                                                os.system("""samtool sort {} > {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+".bam",self.results+map+"/Alignment/"+SM+"_"+LB+"_SE_"+map+"_sorted.bam"))
            print " Mapping using ",map," software is done"
        #Mark and remove the PCR duplicates
    def PCR_dup_mark_remov(self):
             for map in self.mappers:
                   print ("""
                            ==================================================================
                                {} alignment PCR duplication marking and removal begins!
                            ==================================================================
                          """.format(map))
                   for fl in os.listdir(self.results+map+"/Alignment/"):
                                        os.chdir(self.results+map+"/Alignment/")
                                        if fl.endswith(map+"sorted..bam"):
                                            input=os.join.path(self.results+map+"/Alignment/",fl)
                                            output=os.join.path(self.results+map+"/Alignment/",fl.reaplace(map+"_sorted.bam",map+"_dedup.bam"))
                                            os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} OUTPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true""".format(self.picard,input,output))
                                                os.system("""samtools index {}""".format(self.picard,output))
                    return "PCR duplicates removal is done!"
                        
   # Realign around the indels usng GATK the final output is a sorted indexed bam file
    def realignment(self):          
        
        print ("""
                ==================================================================================================
                   Realignment around the indels using GATK: the final output is a sorted indexed Bam file
                ==================================================================================================
            """)
        for map in self.mappers:
            print map
            for fl in os.listdir(self.results+map+"/Alignment/"):
                  os.chdir(self.results+map+"/Alignment/")
                  if fl.endswith(map+"_dedup.bam"):
                        input=os.path.join(self.results+map+"/Alignment/",fl)
                        output1=os.path.join(self.results+map+"/Alignment/",fl.replace(map+"_sorted_dedup.bam",map+".intervals"))
                        output2=os.path.join(self.results+map+"/Alignment/",fl.replace(map+"_sorted_dedup.bam,map+"_realigned.bam"))
                        output3=os.path.join(self.results+map+"/Alignment/",fl.replace(map+"_realigned_.bam",map+"_realig.bam"))
                        os.system("""java -jar {} -T RealignerTargetCreator -nt 8 -R {} -I {} -o {}""".format(self.GATK,self.ref,input,output1))
                        os.system("""java -jar {} -T IndelRealigner -R {} -I {} -targetIntervals {} -o {}""".format(self.GATK,self.ref,input,output1,output2))
                        os.system("""samtools index {} - > {}""".format(output2,output3))
    
        return "realignment around indels is done!"
    
    #Base quality score recalibration
    def base_qual_recal():                   
        
        for map in self.mappers:
             print ("""
                =================================================================
                    {} alignment base quality score recalibration starts now!
                =================================================================
             """.format(map))                     
             for fl in os.listdir(self.results+map+"/Alignment/"):
                  os.chdir(self.results+map+"/Alignment/")
                  if fl.endswith(map+"_realig.bam"):
                        input=os.path.join(self.results+map+"/Alignment/",fl)
                        output1=os.join.path(self.results+map+"/Alignment/",fl.reaplace("_realigned_"+map+".bam",map+"_recal_data.table"))
                        output2=os.join.path(self.results+map+"/Alignment/",fl.reaplace("_realigned_"+map+".bam","_realig_resorted_"+map+".bam"))
                        os.system("""java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R {} -I {} -knownSites {} -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output1))
                        os.system("""java -jar GenomeAnalysisTK.jar -T PrintReads -R {} -BQSR {} -I {} - | samtool sort {} - | samtools index - > {}""".format(self.ref_genome,input,output1,output2))
        
    
       #Check if there is samples from more than one sequencing library if yes merge them and move them to Final_bam folder otherwise use bam file in Alignment directory with the extension final_mapper.bam
    def merge_bam(self):
        for map in self.mappers:
            if not all(len(multi_libraryPE[k]) == 1 for k in multi_libraryPE): #  not all samples have one value (library id): some comes from one library others from different ones
                continue 
                print ("""            
                      ==============================================
                             PE:Libraries merging begins!
                      ============================================== 
                """)
            #samples that do need merging are moved to a new final folder if sample id has one lib value
                for key,value in multi_libraryPE.iteritems():
                    if len(value)==1:#no merging just move bam file to final_Bam folder
                        os.system("cp "+self.results+map+"/Alignment/"+key+"_"+value+"_PE_final_"+map+".bam"+self.results+map+"/Alignment/Final_Bam/"+key+"_"+value+"_PE_final_"+map+".bam ")
                    else: #merge all the bams into final one                  
                        for fl in os.listdir(self.results+map+"/Alignment/"):
                            os.chdir(self.results+map+"/Alignment/")
                            if fl.endswith("_PE_final_"+map+".bam"):      
                                    multi_library=[ ]
                                    to_be_merged=""
                                    for id in SM:
                                        if id in fl:
                                            multi_library.append(self.results+map+"/Alignment/"+fl)
                                            #when i find the id in the file name write the file name multi
                                            for i in multi_library:
                                                    to_be_merged+=i+" "
                                                    output=self.results+map+"/Alignment/Final_Bam/"+id+"_PE_merged_"+map+".bam"
                                                    cmd="""{}samtools merge {} - | {}samtool sort - | {}samtools index - > {}""".format(self.samtools,to_be_merged,self.samtools,self.samtools,output)  
                                                    os.system(cmd)  
            
        #Same processing for samples with SE reads                             
            if not all(len(multi_librarySE[k]) == 1 for k in multi_librarySE):   
                continue 
                print ("""            
                      ==============================================
                             SE:Libraries merging begins!
                      ============================================== 
                """)
                for key,value in multi_librarySE.iteritems():
                    if len(value)==1:
                        os.system("cp "+self.results+map+"/Alignment/"+key+"_"+value+"_SE_final_"+map+".bam"+self.results+map+"/Alignment/Final_Bam/"+key+"_"+value+"_SE_final_"+map+".bam ")
                    else:                   
                        for fl in os.listdir(self.results+map+"/Alignment/"):
                            os.chdir(self.results+map+"/Alignment/")
                            if fl.endswith("_SE_final_"+map+".bam"):      
                                    multi_library=[ ]
                                    to_be_merged=""
                                    for id in SM:
                                      if id in fl:
                                        multi_library.append(self.results+map+"/Alignment/"+fl)
                                        #when i find the id in the file name write the file name multi
                                        for i in multi_library:
                                                    to_be_merged+=i+" "
                                                    output=self.results+map+"/Alignment/Final_Bam/"+id+"_SE_merged_"+map+".bam"
                                                    cmd="""{}samtools merge {} - | {}samtool sort - | {}samtools index - > {}""".format(self.samtools,to_be_merged,self.samtools,self.samtools,output)  
                                                    os.system(cmd)                                          
            else:#when no files needs to be merged exit this function and use 
                exit()           
                                                     
        return " " 

      # Check the mapping quality only samples with >=90% mapped reads and/or genome coverage>=40% are accepted 
    def genome_coverag_mappability_stat(self):
            print ("""
                 ==================================================================
               Genome coverage and reads mappability statistics are starting !
             ==================================================================
                """)            
  
            for map in self.mappers:
                genom_cov={}#sample with genome coverage >40
                read_map={}#sample with mapped reads >90
                if not os.path.isdir(self.results+map+"/statistics/final_stat"):
                    os.makedirs(self.results+map+"/statistics/final_stat")   
                #check if there are merged Bam files or not and call the corresponding folder    
                if os.path.isdir(self.results+map+"/Alignment/Final_Bam"):
                    Final=self.results+map+"/Alignment/Final_Bam"
                else:
                    Final=self.results+map+"/Alignment/"                    
            for fl in os.listdir(Final):
                        os.chdir(Final)      
                    #if fl.endswith("_merged_bwa.bam") or fl.endswith("_final_bwa.bam")      
                        input=Final+"/"+fl
                        for ID in Sample:
                            output1=self.results+map+"/statistics/"+ID+"_bwa_genomecov_bed.txt"
                                        #output2="+map+"_map_output+"/statistics/"+ID+"_bwa_genomecov>40.txt"
                            output2=self.results+map+"/statistics/"+ID+"_bwa_genomecov=0.txt"
                            os.system("""{}bedtools genomecov -ibam -bga {} -g {} > {}""".format(self.bedtools,input,self.ref_genome,output1))
                            cov=os.system("awk '{ print $4 }' "+output1)
                            genom_cov[ID]=cov          
                                        #os.system("""awk 'NF && $4>40' {} > {}""".format(output1,output2))
                            os.system("""awk 'NF && $4<2' {} > {}""".format(output1,output3))
                            output11=self.results+map+"/statistics/"+ID+"_stats_"+map+".txt"
                            os.system("""{}samtools flagstat {} > {}""".format(self.samtools,input,output11))
                            for fl in os.listdir(self.results+map+"/statistics/"): 
                                              if fl.endwith("_stats_"+map+".txt"):          
                                                  for line in open(fl).readlines():
                                                      line=line.strip()
                                                      if "mapped (" in line:
                                                              x=line.split("mapped")[1].split(":")[0].split("(")[1].split("%")[0]
                                                              y=float(x)
                                                              read_map[ID]=y                                                          
                            for (i, j), (k, v) in product(genom_cov.items(), read_map.items()):   
                                                        if i==k and j>=90 or v>=40:
                                                            os.system("""cp {} {}""".format(self.results+map+"/statistics/"+ID+"_stats_"+map+".txt",self.results+map+"/statistics/final_stat"+ID+"_stats_"+map+".txt"))   
            return " "
    
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
            if os.path.isdir(self.results+map+"/Alignment/Final_Bam"):
                Final=self.results+map+"/Alignment/Final_Bam"
            else:
                Final=self.results+map+"/Alignment/"                    
                for fl in os.listdir(Final):
                  os.chdir(Final)           
                  input=Final+"/"+fl              
                  for id in sample:
                                output=self.results+map+"/Joint_variants/"+id+"_GATK_snps_indels.vcf" 
                                os.system("""java -jar {}GenomeAnalysisTK.jar -R {} -T HaplotypeCaller -I {} --dbsnp {} -stand_call_conf 30 -stand_emit_conf 10.0 -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output))


            return "" 
    
    #structural variant identification(big deletion) using DELLY                          
    def structural_variation_calling(self):
        print ("""
                  ================================================
                    structural_variation_callingis starting !
                 ================================================
        """)                    
               
        for map in self.mappers:
            if not os.path.isdir(self.results+map+"/struc_Variants"):
                            os.makedirs(self.results+map+"/struct_Variants")
            if os.path.isdir(self.results+map+"/Alignment/Final_Bam"):
                Final=self.results+map+"/Alignment/Final_Bam"
            else:
                Final=self.results+map+"/Alignment/"                    
                for fl in os.listdir(Final):
                           os.chdir(Final)           
                           input=Final+"/"+fl
                           for id in sample:
                                output=self.results+map+"/struc_Variants"+id+"_bwa_SV"
                                os.system("""{}delly call -t DEL -g {} -o {}.bcf {}""".format(self.DElly,self.ref_genome,output,input))
                                os.system("""{}bcftools call -mvO v -o {}.vcf {}.bcf""".format(self.bcftools,output,output))#convert bcf to vcf
                                
        return "" 
    #SNPs and Indels were annotated as intergenic or genic and for amino acid level changes using ANNOVAR
    def Gene_annotation_annovar(self):
        """Create the gene annotation databases:gtf file and genome assembly fasta
           convert the VCF into annovar format
           """
        print ("""
              ====================================
                   Annotation is starting !
              ====================================
        """)
        if not os.path.isdir(self.TBdb):
            os.makedirs(self.TBdb)              
        #download the genome assembly file from Ensembl
        os.system("wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.toplevel.fa.gz   -P "+self.TBdb) 
        os.system("gunzip "+self.TBdb+"/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.toplevel.fa.gz")
        #convert gene definition file into gene prediction file
        os.system("""chmod +x {}""".format(self.gtfToGenePred))
        for fl in os.listdir(self.TBdb):
            if not os.path.exists(self.TBdb+"/H37Rv_refGene.txt")  and not os.path.exists(self.TBdb+"/H37Rv_refGeneMrna.fa"):                
                os.system("""./{} -genePredExt {} {}""".format(self.gtfToGenePred,self.H37RV_gtf,self.TBdb+"/H37Rv_refGene.txt")) #run a bin file (first check if you have lib2 installed  
        #Generate a transcript fasta file from the genome assembly file
                os.system("""perl {}/retrieve_seq_from_fasta.pl --format refGene --seqfile {} {} --out {}""".format(self.annovar,self.ref_genomv,self.TBdb+"/H37Rv_refGene.txt",self.TBdb+"/H37Rv_refGeneMrna.fa"))              
        #Convert the VCF to annovar input format file
        for map in self.mappers:
            for fl in os.listdir(self.results+map+"/Joint_Variants/"):
                os.chdir(self.results+map+"/Joint_Variants/")
                file_name=os.path.splitext(fl)[0]
                input=fl
                output=self.results+map+"/Joint_Variants/"+file_name+".avinput"
                os.system("""{}convert2annovar.pl -format vcf4 {}} > {}.avinput""".format(self.annovar,input,output))
                os.system("""python -c "import sys;lines=sys.stdin.read();print lines.replace('gi|444893469|emb|AL123456.3|','Chromosome')" < {} > """ .format(output))           
        #run annovar gene annotation():
            if not os.path.isdir(self.results+map+"/Annotation"):
                os.makedirs(self.results+map+"/Annotation")
                for fl in os.listdir(self.results+map+"/Joint_Variants/"):
                    if fl.endswith(".avinput"):
                      os.chdir(self.results+map+"/Joint_Variants/")
                      for id in sample:
                        input=self.results+map+"/Joint_Variants/"+id+".vcf"
                        output=self.results+map+"/Joint_Variants_annovar/"+id+".annovar"
                        os.system("""perl {}table_annovar.pl --vcfinput {} {} --outfile {} --buildver H37Rv --nastring . --protocol refGene --operation g  --vcfinput""".format(self.annovar,input,self.TBdb,output))
        #annotate_variation.pl -filter -dbtype generic -genericdbfile h34rv_generic.txt -build hg19 -out ex1 h37rv.avinput example/        

        return""
        
def main():
        pipeline=BSMTBGan(working_Dir_"Mycobacterium_tuberculosis_h37rv","Mycobacterium_tuberculosis_h37rv.fa","Fastq.csv", mappers=["BWA","Novoalign"])
        try:
            #print pipeline.process_ref_genome()
            #print pipeline.spoligotyping()
            #print pipeline.trimming()
            #print pipeline.Mapping()
            #print pipeline.PCR_dup_mark_remov()
            print pipeline.realignment()
            #print pipeline.base_qual_recal()
            #print pipeline.merge_bam()
            #print pipeline.genome_coverag_mappability_stat()
            #print pipeline.joint_variant_calling()
            #print pipeline.structural_variation_calling()
            #print pipeline.Gene_annotation_annovar()
        except IOError as e:
            print("I/O error: {0}".format(e))    
    
if __name__ == '__main__':main()


