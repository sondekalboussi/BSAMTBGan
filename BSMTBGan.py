"""Please make sure to install java and perl 
converting the gtf file into genepred file must be run on linux OS (no mac/windows version)
self.fastQFileData is a table that contains columns in the following order :file path, ID, PL, SM and LB corresponding identifiers needed to add readgroups to bam files
"""
import os
import sys
import csv
from itertools import product
from datetime import datetime
start_time = datetime.now()
new=open(path to txt file for time recording,"w")
class BSMTBan(object):
    def __init__(self,Novoalign,BWA_mem,BWA_mem_map_output,novo_align_map_output,ref_genome,ref_gen_Dir,Trimmomatic,bwa,novoalign,GATK,fastq_Dir,dbSNP,illumina_adapters,samtools,picard,bedtools,bcftools,Delly,annovar,gtfToGenePred,MTdb)
      #Tool settings and paths
      self.Novoalign=Novoalign
      self.BWA_mem=BWA_mem
      self.BWA_mem_map_output=BWA_mem_map_output
      self.novo_align_map_output=novo_align_map_output
      self.ref_genome=genome#only the reference genome fil
      self.ref_gen_Dir=ref_gen_Dir
      self.trimmomatic=Trimmomatic
      self.bwa=bwa
      self.novoalign=novoalign
      self.GATK=GATk
      self.fastq_Dir=fastq_Dir
      self.dbSNp=dbSNP
      self.illumina_adapters=illumina_adapters
      self.samtools=samtools
      self.picard=picard
      self.bedtools=bedtools
      self.bcftools=bcftools
      self.Delly=Delly
      self.annovar=annovar
      self.fastQFileData=input_table with table 
      self.gtfToGenePred=gtfToGenePred
      self.annotation_db=MTdb
    
    print
	print "_____________________________________________________________________________________________"
	print "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
	print "BSMTBGan - Belgian SouthAfrican MTB genome Annotation pipeline                   "
	print "V1.0, Created by Sondes kalboussi, May 2017, email: sondes777@hotmail.com         "
	print "Copyright Sondes kalboussi 2017                                                           "
	print "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
	print "_____________________________________________________________________________________________"
	print ""
    
    def process_ref_genome(self):
            print ("""
            ==================================================================
			   Reference genomes processing for the annotation is starting"!
			==================================================================
            """)        
        #prepare the ref genome to be used for the mapping and GATK
            for ref_genome in os.listdir(self.ref_gen_Dir):
                os.chdir(self.ref_gen_Dir+"/"+ref_genome+"/FASTA")
                for fl in (self.ref_gen_Dir+"/"+ref_genome+"/FASTA"):
                  if fl.endswith(".fa") or fl.endswith(".fasta"):
                        ID=os.path.splitext(fl)[0]
                        os.system("{}/bwa index {}""".format(self.bwa,fl))
                        os.system("{}/novoindex {}.nix {} "format(self.novoalign,ID,fl))         
                        os.system("{}}/samtools faidx {}""".format(self.samtools,fl))
                        os.system("java -jar {}/picard.jar CreateSequenceDictionary R={}.fasta O={}.dict""".format(self.picard,fl,ID))
            return "Ready to start mapping!"  
          
    def trimming_mapping(self):
        global Sample,Library
        Sample=[]
        Library=[]
        file=csv.reader(open(self.fastQFileData,'rb'))
        next(file, None)
        if self.BWA_mem:
            print ("""
            ========================================================
			 The trimming and BWA mem mapping start now!
			========================================================
            """)
            if not os.path.isdir(self.BWA_mem_map_output+"/Trimming"):
                    os.makedirs(self.BWA_mem_map_output+"/Trimming") 
            if not os.path.isdir(self.BWA_mem_map_output+"/Alignment"):
                os.makedirs(self.BWA_mem_map_output+"/Alignment")    
            for row in file:
                SM=row[3]
                if SM not in Sample:
                  Sample.append(SM)
                ID=row[1]
                PL=row[2]
                if "R1" not in row[0] and "R2" not in row[0]:
                    readSE=row[0]
                    LB=row[4]
                    output_trim=self.BWA_mem_map_output+"/Trimming/"+SM+"_"+LB+"_SE_trimed.fq.gz"
                    os.system("""java -jar {}/trimmomatic.jar SE -threads 8 -phred33 {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(Trimmomatic,readSE,output_trim))
                    readGroup="'@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL'"         
                    output_map=BWA_mem_map_output+"/Alignment/"+SM+"_"+LB"_SE_bwa.bam")
                    os.system("""{}/bwa mem -M -t 16 -R {} {}.fasta {} - | {}/samtools view -Sb - | {}/samtool sort - | {}/samtools index - > {}""".format(self.bwa,readGroup,self.ref_genome,output_trim,self.samtools,self.samtools,self.samtools,output_map))
                if "R1" or "R2" in row[0] and id in row[0] :  
                    read=row[0]
                    lib=row[4]
                    if "R1" in read:
                        read1=read
                        LB1=row[4]
                    if "R2"in read:
                        read2=read
                        LB2=row[4]
                        if Lib1==Lib2:
                          output1=self.BWA_mem_map_output+"/Trimming/"+SM+"R1_"+Lib1+"_PE_paired_trimed.fastq"
                          output2=self.BWA_mem_map_output+"/Trimming/"+SM+"R1_"+Lib1"_PE_unpaired_trimed.fastq"
                          output3=self.BWA_mem_map_output+"/Trimming/"+SM+"R2_"+Lib1"_PE_paired_trimed.fastq"
                          output4=self.BWA_mem_map_output+"/Trimming/"+SM+"R2_"+Lib1"_PE_unpaired_trimed.fastq"
                          readGroup="'@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL'"
                          output_map=self.BWA_mem_map_output+"/Alignment/"+SM+"_"+LB1"_PE_bwa.bam")
                          os.system("""java -jar {}/trimmomatic.jar PE -threads 8 -phred33 {} {} {} {} {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(Trimmomatic,read1,read2,output1,output2,output3,output4))
                          os.system("""{}/bwa mem -M -t 16 -R {} {}.fasta {} {} - | {}/samtools view -Sb - | {}/samtool sort - | {}/samtools index - > {}""".format(self.bwa,readGroup,self.ref_genome,output1,output3,self.samtools,self.samtools,self.samtools,output_map))        
        if self.Novoalign:
            print ("""
            ========================================================
			 The trimming and Novoalign mapping start now!
			========================================================
            """)
            if not os.path.isdir(self.novo_align_map_output+"/Trimming"):
                    os.makedirs(self.novo_align_map_output+"/Trimming")  
            for row in file:
                SM=row[3]
                ID=row[1]
                PL=row[2]
                if "R1" not in row[0] and "R2" not in row[0]:
                    readSE=row[0]
                    LB=row[4]
                    output_trim=self.novo_align_map_output+"/Trimming/"+SM+"_"+LB+"_SE_trimed.fq.gz"
                    os.system("""java -jar {}/trimmomatic.jar SE -threads 8 -phred33 {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(Trimmomatic,readSE,output)) 
                    readGroup="'@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL'"         
                    output_map=self.novo_align_map_output+"/Alignment/"+SM+"_"+LB"_SE_novo.bam")
                    os.system("""{}/bwa mem -M -t 16 -R {} {}.fasta {} - | {}/samtools view -Sb - | {}/samtool sort - | {}/samtools index - > {}""".format(self.novoalign,readGroup,self.ref_genome,output_trim,self.samtools,self.samtools,self.samtools,output_map))
                if "R1" or "R2" in row[0] and id in row[0] :  
                    read=row[0]
                    lib=row[4]
                    if "R1" in read:
                        read1=read
                        LB1=row[4]
                    if "R2"in read:
                        read2=read
                        LB2=row[4]
                        if Lib1==Lib2:
                          output1=self.novo_align_map_output+"/Trimming/"+SM+"R1_"+LB1+"_PE_paired_trimed.fastq"
                          output2=self.novo_align_map_output+"/Trimming/"+SM+"R1_"+LB1"_PE_unpaired_trimed.fastq"
                          output3=self.novo_align_map_output+"/Trimming/"+SM+"R2_"+LB1"_PE_paired_trimed.fastq"
                          output=self.novo_align_map_output+"/Trimming/"+SM+"R2_"+LB1"_PE_unpaired_trimed.fastq"
                          os.system("""java -jar {}/trimmomatic.jar PE -threads 8 -phred33 {} {} {} {} {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(Trimmomatic,read1,read2,output1,output2,output3,output4)) 
                          readGroup="'@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL'"
                          output_map=self.novo_align_map_output+"/Alignment/"+SM+"_"+LB1"_PE_novo.bam"
                          os.system("""{}/bwa mem -M -t 16 -R {} {}.fasta {} {} - | {}/samtools view -Sb - | {}/samtool sort - | {}/samtools index - > {}""".format(self.novoalign,readGroup,self.ref_genome,output1,output3,self.samtools,self.samtools,self.samtools,output_map))
        return "Trimming and mapping using BWA mem and Novolaign are done!"
      
    def realignment(self):          
        #realign around the indels usng GATK the final output is a sorted indexed bam file
          if self.BWA_mem:
              print ("""
            ==================================================================================================
			 BWA_mem: realignment around the indels usong GATK, the final output is a sorted indexed Bam file"
			==================================================================================================
            """)
              for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                  os.chdir(self.BWA_mem_map_output+"/Alignment/")
                  if fl.endswith("_bwa.bam"):
                  		input=os.path.join(self.BWA_mem_map_output+"/Alignment/",fl)
                        output1=os.path.join(self.BWA_mem_map_output+"/Alignment/",fl.reaplace("_bwa.bam","_bwa.intervals"))
                        output2=os.path.join(self.BWA_mem_map_output+"/Alignment/",fl.reaplace("_bwa.bam","_realigned_bwa.bam"))
                        #output3=os.path.join(self.BWA_mem_map_output+"/Alignment/",fl.reaplace("_bwa.bam","_realig_bwa.bam"))
                        os.system("""java -jar {}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {}.fasta -I {} -o {}""".format(self.GATK,self.ref_genome,input,output1)
                        os.system("""java -jar {}/GenomeAnalysisTK.jar -T IndelRealigner -R {}.fasta -I {} -o {} -targetIntervals {} - | {}/samtools index {} - > {}""".format(self.GATK,self.ref_genome,input,input1,self.samtools,output2))
                        #os.system("""{}/samtools index {} - > {}""".format(self.samtools,output2,output3))                  
          if self.Novoalign:
              print ("""
            ==================================================================================================
			Novoalign: realignment around the indels usong GATK, the final output is a sorted indexed Bam file
			==================================================================================================
            """)                    
              for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                  os.chdir(self.novo_align_map_output+"/Alignment/")
                  if fl.endswith("_novo.bam"):                
                  input=os.path.join(self.novo_align_map_output+"/Alignment/",fl)
                  output1=os.path.join(self.novo_align_map_output+"/Alignment/",fl.reaplace("_novo.bam","_novoalg.intervals"))
                  output2=os.path.join(self.novo_align_map_output+"/Alignment/",fl.reaplace("_novo.bam","_realigned_novoalg.bam"))
                  #output3=os.path.join(self.novo_mem_map_output+"/Alignment/",fl.reaplace("_bwa.bam","_realig_bwa.bam"))                
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {}.fasta -I {} -o {}""".format(self.GATK,self.ref_genome,input,output1)
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T IndelRealigner -R {}.fasta -I {}  -targetIntervals {} - | {}/samtools index {} - > {}""".format(self.GATK,self.ref_genome,input,input1,self.samtools,output2))
                  #os.system("""{}/samtools index {} - > {}""".format(self.samtools,output2,output3))   
         return "realignment around indels is done!"  
                            
    def base_qual_recal():                   
        #base quality score recalibration
          if self.BWA_mem:
             print ("""
            =================================================================
			  BWA_mem alignment base quality score recalibration starts now!
			=================================================================
            """)                     
              for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                  os.chdir(self.BWA_mem_map_output+"/Alignment/")
                  if fl.endswith("_realig_bwa.bam"):          
                        input=os.path.join(self.BWA_mem_map_output+"/Alignment/",fl)
                        output1=os.join.path(self.BWA_mem_map_output+"/Alignment/",fl.reaplace("_realig_bwa.bam","bwa_recal_data.table"))
                        output2=os.join.path(self.BWA_map_output+"/Alignment/",fl.reaplace("_realig_bwa.bam","_realig_resorted_bwa.bam"))    
                        os.system("""java -jar {}/GenomeAnalysisTK.jar -T BaseRecalibrator -R {}.fasta -I {} -knownSites {} -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output1))
                        os.system("""java -jar {}/GenomeAnalysisTK.jar -T PrintReads -R {}.fasta -BQSR {} -I {} - | {}/samtool sort {} - | {}/samtools index - > {}""".format(self.GATK,self.ref_genome,input,output1,output2))
         if self.Novoalign:
              print ("""
            ==================================================================
			  Novoalign alignment base quality score recalibration starts now!
			==================================================================
            """)                             
              for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                  os.chdir(self.novo_align_map_output+"/Alignment/")
                  if fl.endswith("_realig_novoalg.bam"):                            
                        input=os.join.path(self.novo_align_map_output+"/Alignment/",fl)
		               	output1=os.join.path(self.novo_align_map_output+"/Alignment/",fl.replace("_realig_novoalg.bam","_novoalg_recal_data.table"))
                        output2=os.join.path(self.novo_align_map_output+"/Alignment/",fl.reaplace("_realig_novo.bam","_realig_resorted_novo.bam"))    
                        os.system("""java -jar {}/GenomeAnalysisTK.jar -T BaseRecalibrator -R {}.fasta -I {} -knownSites {} -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output1))
                        os.system("""java -jar {}/GenomeAnalysisTK.jar -T PrintReads -R {}.fasta -BQSR {} -I {} - | {}/samtool sort {} - | {}/samtools index - > {}""".format(self.GATK,self.ref_genome,input,output1,output2))
         return "Base quality score recalibration is done!""  
                            
    def PCR_dup_mark_remov(self):            
        #mark the duplicates
            if self.BWA_mem:
                 print ("""
            ==================================================================
			  BWA_mem alignment PCR duplication marking and removal begins!
			==================================================================
            """)                      
                for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):         
                    os.chdir(self.BWA_mem_map_output+"/Alignment/")
                    if fl.endswith("_realig_resorted_bwa.bam"):        
                          input=os.join.path(self.BWA_mem_map_output+"/Alignment/",fl)
                          output=os.join.path(self.BWA_map_output+"/Alignment/",fl.reaplace("_realig_resorted_bwa.bam","_final_bwa.bam"))
                            
                          #os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} OUTPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true""".format(self.picard,input,output))
                          os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true - | {}/samtools index {} - > {} """.format(self.picard,input,self.samtools,output))  
                          #os.system("""{}/samtools index {} - > {}""".format(self.samtools,output))
            if self.Novoalign:
                print ("""
            ==================================================================
			  Novoalign alignment PCR duplication marking and removal begins!
			==================================================================
            """)              
                for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                    os.chdir(self.novo_align_map_output+"/Alignment/")
                    if fl.endswith("_realig_resorted_novo.bam"):          
                            input=os.join.path(self.novo_align_map_output+"/Alignment/",fl)
                            output=os.join.path(self.novo_align_map_output+"/Alignment/",fl.reaplace("_realig_resorted_novo.bam","_final_novo.bam"))
                            #os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} OUTPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true""".format(self.picard,input,output))
                            os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true - | {}/samtools index {} - > {} """.format(self.picard,input,self.samtools,output)) 
                            #os.system("""{}/samtools index {}""".format(self.samtools,output))
            return "PCR duplicates removal is done!"
                            
    def merge_bam(self):
        # merge bam files of same sample from different library and sort, index the merged bam                   
        multi_libraryPE=[ ]
        multi_librarySE=[ ] 
        for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                    os.chdir(self.BWA_mem_map_output+"/Alignment/")
                    if fl.endswith("_PE_final_bwa.bam"):      
                                IDPE=fl.split("_PE_final_bwa.bam")[0]
                                multi_libraryPE.append(IDPE)
                    if fl.endswith("_SE_final_bwa.bam"):
                                IDSE=fl.split("_PE_final_bwa.bam")[0]
                                multi_librarySE.append(IDSE)
         for element in (multi_libraryPE):
             if sample.count(element)==1:
                print ("""            
                    =====================================================================
                     No Bam files to be merged!
                    =====================================================================  
                    """)            
                break
                            
             else:
                  if self.BWA_mem:        
                       print ("""            
                      =====================================================================
                       BWA_mem BAM paired end files merging begins!
                      =====================================================================  
                      """)
                       for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                              os.chdir(self.BWA_mem_map_output+"/Alignment/")       
                              if fl.endswith("_PE_final_bwa.bam"):      
                                  multi_library=[ ]
                                  to_be_merged=""
                                  for id in SM:
                                    if id in fl:
                                        multi_library.append(self.BWA_mem_map_output+"/Alignment/"+fl)
                                    #when i find the id in the file name write the file name multi
                                        for i in multi_library:
                                            to_be_merged+=i+" "
                                            output=self.BWA_map_output+"/Alignment/"+id+"_PE_merged_bwa.bam"
                                            cmd="""{}/samtools merge {} - | {}/samtool sort - | {}/samtools index - > {}""".format(self.samtools,to_be_merged,self.samtools,self.samtools,output)  
                                            os.system(cmd)  
                            
                    if self.Novoalign: 
                          print ("""            
                          =====================================================================
                           Novoalign BAM paired end files merging begins!
                          =====================================================================  
                          """)       
 
                           for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                                            os.chdir(self.novo_align_map_output+"/Alignment/")
                                            if fl.endswith("_PE_final_bwa.bam"):      
                                                  multi_library=[ ]
                                                  to_be_merged=""
                                                  for id in SM:
                                                    if id in fl:
                                                        multi_library.append(self.BWA_mem_map_output+"/Alignment/"+fl)
                                                    #when i find the id in the file name write the file name multi
                                                        for i in multi_library:
                                                            to_be_merged+=i+" "
                                                            output=self.BWA_map_output+"/Merged_alignment/"+id+"_PE_merged_bwa.bam"
                                                            cmd="""{}/samtools merge {} - | {}/samtool sort - | {}/samtools index - > {}""".format(self.samtools,to_be_merged,self.samtools,self.samtools,output)  
                                                            os.system(cmd)         
                            
         for element in (multi_librarySE):
             if sample.count(element)==1:
                print ("""            
                    =====================================================================
                     No Bam files to be merged!
                    =====================================================================  
                    """) 
                break            
            else:
                if self.BWA_mem:        
                       print ("""            
                      =====================================================================
                       BWA_mem BAM single end files merging begins!
                      =====================================================================  
                      """)            
                      for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                                  os.chdir(self.BWA_mem_map_output+"/Alignment/")                
                                  if fl.endswith("_SE_final_bwa.bam"):      
                                      multi_library=[ ]
                                      to_be_merged=""
                                      for id in SM:
                                        if id in fl:
                                            multi_library.append(self.BWA_mem_map_output+"/Alignment/"+fl)
                                        #when i find the id in the file name write the file name multi
                                            for i in multi_library:
                                                to_be_merged+=i+" "
                                                output=self.BWA_map_output+"/Merged_alignment/"+id+"_SE_merged_bwa.bam"
                                                cmd="""{}/samtools merge {} - | {}/samtool sort - | {}/samtools index - > {}""".format(self.samtools,to_be_merged,self.samtools,self.samtools,output)  
                                                os.system(cmd)                             
                if self.Novoalign:
                       print ("""            
                      =====================================================================
                       Novoalign BAM single end files merging begins!
                      =====================================================================  
                      """)      
                      for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                                            os.chdir(self.novo_align_map_output+"/Alignment/")
                  					  if fl.endswith("_SE_final_bwa.bam"):      
                                                  multi_library=[ ]
                                                  to_be_merged=""
                                                  for id in SM:
                                                    if id in fl:
                                                        multi_library.append(self.BWA_mem_map_output+"/Alignment/"+fl)
                                                    #when i find the id in the file name write the file name multi
                                                        for i in multi_library:
                                                            to_be_merged+=i+" "
                                                            output=self.BWA_map_output+"/Alignment/"+id+"_SE_merged_bwa.bam"
                                                            cmd="""{}/samtools merge {} - | {}/samtool sort - | {}/samtools index - > {}""".format(self.samtools,to_be_merged,self.samtools,self.samtools,output)  
                                                            os.system(cmd)                                   
            return "merging bam is done!" 

    def genome_coverag_mappability_stat(self):
  
            if self.BWA_mem:
                genom_cov={}#sample with genome coverage >40
                read_map={}#sample with mapped reads >90
                if not os.path.isdir(self.BWA_mem_map_output+"/statistics/final_stat"):
                    os.makedirs(self.BWA_mem_map_output+"/statistics/final_stat")   
                for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                    os.chdir(self.BWA_mem_map_output+"/Alignment/")
                    		if fl.endswith("_merged_bwa.bam") or fl.endswith("_final_bwa.bam")      
                                  input=self.BWA_mem_map_output+"/Alignment/"+fl
                                  for ID in Sample:
                                        output1=self.BWA_mem_map_output+"/statistics/"+ID+"_bwa_genomecov_bed.txt"
                                        #output2=BWA_mem_map_output+"/statistics/"+ID+"_bwa_genomecov>40.txt"
                                        output2=self.BWA_mem_map_output+"/statistics/"+ID+"_bwa_genomecov=0.txt"
                                        os.system("""{}/bedtools genomecov -ibam -bga {} -g {}.fasta > {}""".format(self.bedtools,input,self.ref_genome,output1)
                                        cov=os.system("awk '{ print $4 }' "+output1)
                                        genom_cov[ID]=cov          
                                        #os.system("""awk 'NF && $4>40' {} > {}""".format(output1,output2))
                                        os.system("""awk 'NF && $4<2' {} > {}""".format(output1,output3))
                                        output11=self.BWA_mem_map_output+"/statistics/"+ID+"_stats_bwa.txt"
                                        os.system("""{}/samtools flagstat {} > {}""".format(self.samtools,input,output11))
                                        for fl in os.listdir(self.BWA_mem_map_output+"/statistics/"): 
                                              if fl.endwith("_stats_bwa.txt"):          
                                                  for line in open(fl).readlines():
                                                      line=line.strip()
                                                      if "mapped (" in line:
                                                              x=line.split("mapped")[1].split(":")[0].split("(")[1].split("%")[0]
                                                              y=float(x)
                                                              read_map[ID]=y                                                          
                                          for (i, j), (k, v) in product(genom_cov.items(), read_map.items()):   
                                                        if i==k and j>=90 or v>=40:
                                                            os.system("""cp {} {}""".format(self.BWA_mem_map_output+"/statistics/"+ID+"_stats_bwa.txt",self.BWA_mem_map_output+"/statistics/final_stat"+ID+"_stats_bwa.txt"))   
            if self.Novoalign:
                genom_cov={}
                read_map={}
                if not os.path.isdir(self.novo_align_map_output+"/statistics/final_stat"):
                    os.makedirs(self.novo_align_map_output+"/statistics/final_stat")   
                for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                    os.chdir(self.novo_align_map_output+"/Alignment/")
					if fl.endswith("_merged_bwa.bam") or fl.endswith("_final_bwa.bam")      
                                  input=self.BWA_mem_map_output+"/Alignment/"+fl
                                  for ID in Sample:
                                          input=self.novo_align_map_output+"/Alignment/"+fl                    
                                          output1=self.novo_align_map_output+"/statistics/"+ID+"_novo_genomecov_bed.txt"
                                          #output2=novo_align_map_output+"/statistics/"+ID+"_novo_genomecov>40.txt"#genome coverage 40 or more
                                          output2=self.novo_align_map_output+"/statistics/"+ID+"_novo_genomecov=0.txt"#genome coverage equal to zero
                                          os.system("""{}/bedtools genomecov -ibam -bga {} -g {}.fasta > {}""".format(self.bedtools,input,self.ref_genome,output1)
                                          cov=os.system("awk '{ print $4>=40 }' "+output1)
                                          genom_cov[ID]=cov          
                                          #os.system("""awk 'NF && $4>=40' {} > {}""".format(output1,output2))
                                          os.system("""awk 'NF && $4<2' {} > {}""".format(output1,output2))
                                          output11=self.novo_align_output+"/statistics/"+ID+"_stats_novo.txt"
                                          os.system("""{}/samtools flagstat {} > {}""".format(self.samtools,input,output11))
                                          for fl in os.listdir(self.novo_align_mem_map_output+"/statistics"):
                                             if fl.endwith("_stats_novo.txt"):           
                                              for line in open(fl).readlines():
                                                  line=line.strip()
                                                  if "mapped (" in line:
                                                          x=line.split("mapped")[1].split(":")[0].split("(")[1].split("%")[0]
                                                          y=float(x)
                                                          read_map[ID]=y
                                          for (i, j), (k, v) in product(genom_cov.items(), read_map.items()):   
                                                        if i==k and j>=90 or v>=40:                    
                                                            os.system("""cp {} {}""".format(self.novo_align_map_output+"/statistics/"+ID+"_stats_novo.txt",self.novo_align_map_output+"/statistics/final_stat/"+ID+"_stats_novo.txt"))                          
            print self.BWA_mem_map_output+"/statistics/final_stat"
            print self.novo_align_map_output+"/statistics/final_stat"
            return ""

    def joint_variant_calling(self):
               # joint variant call SNP and indels using GATK
                if self.BWA_mem:
                    if not os.path.isdir(self.BWA_mem_map_output+"/Joint_Variants"):
                            os.makedirs(self.BWA_mem_map_output+"/Joint_Variants")
                    for fl in os.listdir(self.BWA_mem_map_output+"/Merged_alignment/"):
                            os.chdir(self.BWA_mem_map_output+"/Merged_alignment/")
                            for id in sample:
                                input=self.BWA_mem_map_output+"/Fianl_alignment/"+id+"_merged_bwa.bam"
                                output=self.BWA_map_output+"Variants/"+id+"_GATK_snps_indels.vcf" 
                                os.system("""java -jar {}/GenomeAnalysisTK.jar -R {}.fasta -T HaplotypeCaller -I {} --dbsnp {} -stand_call_conf 30 -stand_emit_conf 10.0 -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output,))

                if self.Novoalign:
                    if not os.path.isdir(self.novo_align_map_output+"/Joint_Variants"):
                            os.makedirs(self.novo_align_map_output+"/Joint_Variants")
                    for fl in os.listdir(self.novo_align_map_output+"/Merged_alignment"):
                        os.chdir(self.novo_align_map_output+"/Merged_alignment")
                        for id in sample:
                            input=self.novo_align_map_output+"/Merged_alignment/"+id+"_merged_novo.bam"
                            output=self.novo_align_output+"/reads_stat/"+id+"_GATK_snps_indels.vcf"
                            os.system("""java -jar {}/GenomeAnalysisTK.jar -R {}.fasta -T HaplotypeCaller -I {} --dbsnp {} -stand_call_conf 30 -stand_emit_conf 10.0 -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output,))
                
                print self.BWA_mem_map_output+"/Variants"
                print self.novo_align_map_output+"/Variants"
                return "" 
        
    def structural_variation_calling(self):
                #structural variant (big deletion) using DELLY
                if self.BWA_mem:
                    if not os.path.isdir(self.BWA_mem_map_output+"/struc_Variants"):
                            os.makedirs(self.BWA_mem_map_output+"/struct_Variants")
                    for fl in os.listdir(self.BWA_mem_map_output+"/Merged_alignment/"):
                            os.chdir(self.BWA_mem_map_output+"/Merged_alignment/")
                            for id in sample:
                                input=self.BWA_mem_map_output+"/Fianl_alignment/"+id+"_merged_bwa.bam"#sorted indexed
                                output=self.BWA_mem_map_output+"/struc_Variants"+id+"_bwa_SV"
                                os.system("""{}/delly call -t DEL -g {} -o {}.bcf {}""".format(self.DElly,self.ref_genome,output,input)
                                os.system("""{}/bcftools call -mvO v -o {}.vcf {}.bcf""".format(self.bcftools,output,output))#convert bcf to vcf
                                
                if self.Novoalign:
                    if not os.path.isdir(self.novo_align_map_output+"/struc_Variants"):
                            os.makedirs(self.novo_align_map_output+"/struct_Variants")
                    for fl in os.listdir(self.novo_align_map_output+"/Merged_alignment/"):
                            os.chdir(self.novo_align_map_output+"/Merged_alignment/")
                            for id in sample:
                                input=self.novo_align_map_output+"/Fianl_alignment/"+id+"_merged_bwa.bam"#sorted indexed
                                output=self.novo_align_map_output+"/struc_Variants"+id+"_novo_SV"
                                os.system("""{}/delly call -t DEL -g {} -o {}.bcf {}""".format(self.DElly,self.ref_genome,output,input)
                                os.system("""{}/bcftools call -mvO v -o {}.vcf {}.bcf""".format(self.bcftools,output,output))#convert bcf to vcf                      
                print self.BWA_mem_map_output+"/struc_Variants"
                print self.novo_align_map_output+"/struc_Variants"
                return "" 
        
    def mapping_functional_annotation_annovar(self):
             #prepare annotation_files: convert gtf into genepred file
                for fl in os.listdir(self.annotation_db):
                     if not os.path.exists(self.annotation_db+"/H37RV_refGene.txt")                   
                        os.system("""{} -genePredExt {} """.format(self.gtfToGenePred,self.H37RV.gtf,self.annotation_db+"/H37RV_refGene.txt")                             
             #run annovar annotation():
                if self.BWA_mem:
                    if not os.path.isdir(self.BWA_mem_map_output+"/Annotation"):
                            os.makedirs(self.BWA_mem_map_output+"/Annotation")
                    for fl in os.listdir(self.BWA_mem_map_output+"/Joint_Variants/"):
                            os.chdir(self.BWA_mem_map_output+"/Joint_Variants/")
                            for id in sample:
                                input=self.BWA_mem_map_output+"/Joint_Variants/"+id+".vcf"
                                output=self.BWA_mem_map_output+"/Joint_Variants_annovar/"+id+".annovar"
                                os.system("""perl {}/table_annovar.pl --vcfinput {} {} --outfile {} --buildver H37RV --nastring . --protocol refGene,dbsnp --operation g,r,f  --vcfinput""".format(self.annovar,input,TBdb,output)
            
                if self.Novoalign:
                    if not os.path.isdir(self.novo_align_map_output+"/Annotation"):
                            os.makedirs(self.novo_align_map_output+"/Annotation")
                    for fl in os.listdir(self.novo_align_map_output+"/Joint_Variants/"):
                            os.chdir(self.novo_align_map_output+"/Joint_Variants/")
                            for id in sample:
                                input=self.novo_align_map_output+"/Joint_Variants/"+id+".vcf"
                                output=self.novo_align_map_output+"/Annotation/"+id+".annovar"
                                os.system("""perl {}/table_annovar.pl --vcfinput {} {} --outfile {} --buildver H37RV --nastring . --protocol refGene,dbsnp --operation g,r,f  --vcfinput""".format(self.annovar,input,TBdb,output)          
                                
                return""
                
if __name__ == '__main__':
    
    pipeline=mapping_to_ref_genome_functional_annotation(Novoalign,BWA_mem,BWA_mem_map_output,novo_align_map_output,ref_genome,ref_gen_Dir,Trimmomatic,bwa,novoalign,GATK,fastq_Dir,dbSNP,illumina_adapters,samtools,picard,bedtools,bcftools,Delly,annovar,gtfToGenePred,MTdb)                                      
    print pipeline.process_ref_genome()
    print pipeline.run_mapping()
    print pipeline.realignment
    print pipeline.base_qual_recal()
    print pipeline.PCR_dup_remov()
    print pipeline.merge_bam()
    print pipeline.genome_coverag_mappability_stat()
    print pipeline.joint_variant_calling()
    print pipeline.structural_variation_calling()                                      
    print pipeline.mapping_functional_annotation_annovar()

end_time = datetime.now()
new.write('Pipeline execution time: {}'.format(end_time - start_time))
new.close()

