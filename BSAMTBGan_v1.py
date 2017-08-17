"""Please make sure to install java and perl 
converting the gtf file into genepred file must be run on linux OS (no mac/windows version)
self.fastQFileData is a table that contains columns in the following order :file path, ID, PL, SM and LB corresponding identifiers needed to add readgroups to bam files
install libpng12-0 :sudo apt-get install 

"""
import os
import sys
import csv
from itertools import product
from datetime import datetime
start_time = datetime.now()
new=open(path to txt file for time recording,"w")
class BSMTBGan(object):
    def __init__(global_Dir,user_ref,fastq_Dir,File_table,mappers):
      
      #Settings parameters and paths
          self.global_Dir=global_Dir# path to home directory in the cluster or local machine
          self.results=self.global_Dir+"/Results/"#results folder
          self.tools=Tools#Tools binaries in the cluster
          self.ref_gen_Dir=self.global_Dir+"/Software/USAP/Reference/"#H37rv reference
          self.ref_genome=self.ref_gen_Dir+user_ref+"/Fasta/"#fasta file of reference genome
          self.novoalign=self.tools+"/novocraft/3.02.07/bin/"# path to the novocraft software
          self.BWA=self.tools+"/0.7.12/bin/"#path to bwa mem software
          self.trimmomatic="/home/opt/software/Trimmomatic-0.36/"
          self.GATK=self.tools+"/gatk/GATK_3.6.0/"
          self.samtools=self.tools+"/samtools/1.3.1/"
          self.picard=self.tools+"/picard/picard-tools-2.3.0/"
          self.bedtools="/opt/software/BEDTools/2.16.2/"
          self.bcftools=self.tools+"bcftools/bcftools/1.4/"
          self.Delly=self.globalDir+"/Software/BSAMTBan/Tools_new/delly/src/"
          self.annovar=self.tools+"/ANNOVAR/2017-06-01/"
          self.fastq_Dir=self.global_Dir+"/"+fastq_Dir
          self.dbSNp=self.globalDir+"/Software/USAP/Reference/"+userRef+"/dbSNP/"
          self.illumina_adapters=self.globalDir+"/Software/USAP/Tools/illumina_adapters.fna.fasta"
          self.fastQFileData=File_table 
          self.gtfToGenePred=self.globalDir+"/Software/BSAMTBan/Tools_new/dbMTB/"#bin file
          self.H37RV_gtf=#gtf file
          self.TBdb=self.globalDir+"/Software/BSAMTBan/Tools_new/dbMTB/"
          self.spotyping=self.globalDir+"/Software/BSAMTBan/Tools/SpoTyping-2.1/SpoTyping-v2.1-commandLine/"
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
    
    #Process the reference genome for the later use for mapping and GATK
   def process_ref_genome(self):
            print ("""
                ============================================================================
		        Reference genomes processing for the annotation is starting"!
	        ============================================================================
            """)        
      
           for ref_genome in os.listdir(self.ref_gen_Dir):
                os.chdir(self.ref_genome)
               for fl in (self.ref_genome):
                  if fl.endswith(".fa") or fl.endswith(".fasta"):
                        ID=os.path.splitext(fl)[0]
                        os.system("{}bwa index {}""".format(self.BWA,fl))
                        os.system("{}novoindex {}.nix {} "format(self.Novoalign,ID,fl))         
                       os.system("{}samtools faidx {}""".format(self.samtools,fl))
                       os.system("java -jar {}picard.jar CreateSequenceDictionary R={}.fasta O={}.dict""".format(self.picard,fl,ID))
 
	   return " "  

    # Strain identification using Spotyping
    def spoligotyping(self):
        global Sample
        Sample=[]
        #use Spotyping software to classify Mtb strains  
        file=csv.reader(open(self.fastQFileData,'rb'))
        next(file, None)
        print ("""
                  ========================================================
		               Strain identification starts now!
	          ========================================================
        """)
            if not os.path.isdir(self.results+"/Spoligotyping"):
                    os.makedirs(self.results+"/Spoligotyping") 
            for row in file:
                SM=row[3]
                if SM not in Sample:
                    Sample.append(SM)
                if "R1" not in row[0] and "R2" not in row[0]:
                    readSE=row[0]
                    output_spolig_SE=self.results+"/Spoligotyping/"+SM+"_"+LB+"_SE_spotyping.txt"
                    os.system("""python2.7 {}SpoTyping.py {} -o {}""".format(self_spotyping,readSE,output_spolig_SE))
                if "R1" or "R2" in row[0] and SM in row[0] :  
                    read=row[0]
                    if "R1" in read:
                        read1=read
                    if "R2"in read: 
                        read2=read 
                        output_spolig_PE=self.results+"/Spoligotyping/"+SM+"_"+LB+"PE_spotyping.txt"
                        os.system("""python2.7 {}SpoTyping.py {} {} -o {}""".format(self_spotyping,read1,read2,output_spolig_PE))
        return ""   
    
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
                if "R1" not in row[0] and "R2" not in row[0]:
                    readSE=row[0]
                    LB=row[4]
                    output_trim=self.results+"Trimming/"+SM+"_"+LB+"_SE_trimed.fq.gz"
                    os.system("""java -jar {}trimmomatic.jar SE -threads 8 -phred33 {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(self.trimmomatic,readSE,output_trim))
                    readGroup="'@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL'"
                if "R1" or "R2" in row[0] and SM in row[0] :  
                    read=row[0]
                    lib=row[4]
                    if "R1" in read:
                        read1=read
                        LB1=row[4]
                    if "R2"in read:
                        read2=read
                        LB2=row[4]
                        if Lib1==Lib2:
                            output1=self.results+"Trimming/"+SM+"R1_"+Lib1+"_PE_paired_trimed.fq"
                            output2=self.results+"Trimming/"+SM+"R1_"+Lib1"_PE_unpaired_trimed.fq"
                            output3=self.results+"Trimming/"+SM+"R2_"+Lib1"_PE_paired_trimed.fq"
                            output4=self.results+"Trimming/"+SM+"R2_"+Lib1"_PE_unpaired_trimed.fq"
                            readGroup="'@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL'"
                            os.system("""java -jar {}trimmomatic.jar PE -threads 8 -phred33 {} {} {} {} {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(Trimmomatic,read1,read2,output1,output2,output3,output4))
        
	# Create a dictionary where the keys are samples ID and value is the reads ID when they are coming from more than sequencing library so later you can merge all the output for the same sample 
        for fl in os.listdir(self.results+"/Trimming/"):
                        os.chdir(self.results+"/Trimming/")
                        if fl.endswith("_PE_paired_trimed.fq"): 
                                    IDPE=f.split("_PE_paired_trimed.fq")[0]
                                    sm=IDPE.split("_")[0]
                                    lib=IDPE.split("_")[1]
                                    if sm not in multi_libraryPE.keys():
                                        multi_libraryPE[sm]=[lib]
                                    else:
                                        multi_libraryPE[sm].append(lib)
                        if fl.endswith("_SE_final_bwa.bam"):
                                    IDSE=f.split("_SE_final_bwa.bam")[0]
                                    sm=IDSE.split("_")[0]
                                    lib=IDSE.split("_")[1]
                                    if sm not in multi_librarySE.keys():
                                        multi_librarySE[sm]=[lib]
                                    else:
                                        multi_librarySE[sm].append(lib)                  
        return "Trimming is done!"
   
   # Mapping reads using BWA-mem and Novoalign, mappers is a list of software that we want to use 
   def Mapping(self):
	    for map in mappers:# mappers is a list of mapping software name (BWA-mem and novoalign in our study)
		    if not os.path.isdir(self.results+map+"/Alignment"): 
			    os.makedirs(self.results++map+"/Alignment")
		    for fl in os.listdir(self.results+"Trimming/"):
			    os.chdir(self.results+"Trimming/")
                input=os.join.path(self.results+"Trimming/",fl)
		#Mapping the SE reads
                if "SE" in input:
                    output_map=self.results+"/"+map+"/"+SM+"_"+LB"_SE_"+map+".bam")
                    if map=="BWA":
                        os.system("""{}bwa mem -M -t 16 -R {} {}.fasta {} - | {}samtools view -Sb - | {}samtool sort - | {}samtools index - > {}""".format(self.BWA,readGroup,user_ref,output_trim,self.samtools,self.samtools,self.samtools,output_map_bwa))
                    if map=="Novoalign":
                        for fl in os.listdir(self.ref_gen_Dir+user_ref++"/Fasta/"):
                                        os.chdir(self.ref_gen_Dir+user_ref+"/Fasta/"):
                                        if fl.endswith(".nix"):  
                                            os.system("""{}novoalign -d {} -f {} -o SAM | {}samtools view -Sb - | {}samtool sort - | {}samtools index - > {}""".format(self.novoalign,fl,readSE,self.samtools,self.samtools,self.samtools,output_map_novo))
                #Mapping the PE reads
                if "PE" in input:
                    if "R1" in input:
                        R1=input
                    if "R2" in input:
                        R2=input    
                        output_map=self.results+"/"+map+"/"+SM+"_"+LB"_PE_"+map+".bam")
                        if map=="BWA":
                            os.system("""{}bwa mem -M -t 16 -R {} {}.fasta {} {} - | {}samtools view -Sb - | {}samtool sort - | {}samtools index - > {}""".format(self.bwa,readGroup,user_ref,R1,R2,self.samtools,self.samtools,self.samtools,output_map))        
                        if self.Novoalign:
                            for fl in os.listdir(self.ref_gen_Dir+user_ref++"/Fasta/"):
                                os.chdir(self.ref_gen_Dir+user_ref++"/Fasta/"):
                                if fl.endswith(".nix"):    
                                    os.system("""{}novoalign -d {} -f {} {} -o SAM | {}samtools view -Sb - | {}samtool sort - | {}samtools index - > {}""".format(self.novoalign,fl,R1,R2,self.samtools,self.samtools,self.samtools,output_map_novo))
                 print " Mapping using ",map," software is done"
        return " "
   
   # Realign around the indels usng GATK the final output is a sorted indexed bam file
   def realignment(self):          
        for map in mappers:  
            print ("""
                ==================================================================================================
	          {}:realignment around the indels using GATK, the final output is a sorted indexed Bam file"
	        ==================================================================================================
            """.format(map))
            for fl in os.listdir(self.results+"/"+map+"/Alignment/"):
                  os.chdir(self.results+"/"+map+"/Alignment/")
                  if fl.endswith(map+".bam"):
                  		input=os.path.join(self.results+"/"+map+"/Alignment/",fl)
                        output1=os.path.join(self.results+"/"+map+"/Alignment/",fl.reaplace(map+".bam",map+".intervals"))
                        output2=os.path.join(self.results+"/"+map+"/Alignment/",fl.reaplace(map+".bam","_realigned_"+map+".bam"))
                        #output3=os.path.join(self.""+map+"_map_output+"/Alignment/",fl.reaplace("_bwa.bam","_realig_bwa.bam"))
                        os.system("""java -jar {}GenomeAnalysisTK.jar -T RealignerTargetCreator -R {}.fasta -I {} -o {}""".format(self.GATK,user_ref,input,output1)
                        os.system("""java -jar {}GenomeAnalysisTK.jar -T IndelRealigner -R {}.fasta -I {} -o {} -targetIntervals {} - | {}samtools index {} - > {}""".format(self.GATK,user_ref,input,input1,self.samtools,output2))
                        #os.system("""{}/samtools index {} - > {}""".format(self.samtools,output2,output3))                  
    
        return "realignment around indels is done!"  
    
    #Base quality score recalibration                        
    def base_qual_recal():                   
        
        for map in mappers: 
             print ("""
                =================================================================
	           {} alignment base quality score recalibration starts now!
	        =================================================================
             """.format(map))                     
            for fl in os.listdir(self.results+"/"+map+"/Alignment/"):
                  os.chdir(self.results+"/"+map+"/Alignment/")
                  if fl.endswith("_realigned_"+map+".bam"):          
                        input=os.path.join(self.results+"/"+map+"/Alignment/",fl)
                        output1=os.join.path(self.results+"/"+map+"/Alignment/",fl.reaplace("_realigned_"+map+".bam",map+"_recal_data.table"))
                        output2=os.join.path(self.results+"/"+map+"/Alignment/",fl.reaplace("_realigned_"+map+".bam","_realig_resorted_"+map+".bam"))    
                        os.system("""java -jar {}/GenomeAnalysisTK.jar -T BaseRecalibrator -R {}.fasta -I {} -knownSites {} -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output1))
                        os.system("""java -jar {}/GenomeAnalysisTK.jar -T PrintReads -R {}.fasta -BQSR {} -I {} - | {}/samtool sort {} - | {}/samtools index - > {}""".format(self.GATK,self.ref_genome,input,output1,output2))
        return "Base quality score recalibration is done!""  
    
    #Mark and remove the PCR duplicates				  
    def PCR_dup_mark_remov(self):            
        for map in mappers: 
                print ("""
                        ==================================================================
	                   {} alignment PCR duplication marking and removal begins!
	        	==================================================================
                """.format(map))                      
                for fl in os.listdir(self.results+"/"+map+"/Alignment/"):         
                    os.chdir(self.results+"/"+map+"/Alignment/")
                    if fl.endswith("_realig_resorted_"+map+".bam"):        
                          input=os.join.path(self.results+"/"+map+"/Alignment/",fl)
                          output=os.join.path(self.results+"/"+map+"/Alignment/",fl.reaplace("_realig_resorted_"+map+".bam","_final_"+map+".bam"))       
                          #os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} OUTPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true""".format(self.picard,input,output))
                          os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true - | {}/samtools index {} - > {} """.format(self.picard,input,self.samtools,output))  
                          #os.system("""{}/samtools index {} - > {}""".format(self.samtools,output))
            
        return "PCR duplicates removal is done!"
    
    #Check if there is samples from more than one sequencing library if yes merge them and move them to Final_bam folder otherwise use bam file in Alignment directory with the extension final_mapper.bam                       				  
    def merge_bam(self):
        for map in mappers:     
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
                        os.system("cp "+self.results+"/"+map+"/Alignment/"+key+"_"+value+"_PE_final_"+map+".bam"+self.results+"/"+map+"/Alignment/Final_Bam/"+key+"_"+value"_PE_final_"+map+".bam ") 
                    else: #merge all the bams into final one                  
                        for fl in os.listdir(self.results+"/"+map+"/Alignment/"):
                            os.chdir(self.results+"/"+map+"/Alignment/")       
                                if fl.endswith("_PE_final_"+map+".bam"):      
                                    multi_library=[ ]
                                    to_be_merged=""
                                    for id in SM:
                                    if id in fl:
                                        multi_library.append(self.results+"/"+map+"/Alignment/"+fl)
                                            #when i find the id in the file name write the file name multi
                                        for i in multi_library:
                                                    to_be_merged+=i+" "
                                                    output=self.results+"/"+map+"/Alignment/Final_Bam/"+id+"_PE_merged_"+map+".bam"
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
                        os.system("cp "+self.results+"/"+map+"/Alignment/"+key+"_"+value+"_SE_final_"+map+".bam"+self.results+"/"+map+"/Alignment/Final_Bam/"+key+"_"+value"_SE_final_"+map+".bam ") 
                    else:                   
                        for fl in os.listdir(self.results+"/"+map+"/Alignment/"):
                            os.chdir(self.results+"/"+map+"/Alignment/")       
                                if fl.endswith("_SE_final_"+map+".bam"):      
                                    multi_library=[ ]
                                    to_be_merged=""
                                    for id in SM:
                                    if id in fl:
                                        multi_library.append(self.results+"/"+map+"/Alignment/"+fl)
                                        #when i find the id in the file name write the file name multi
                                        for i in multi_library:
                                                    to_be_merged+=i+" "
                                                    output=self.results+"/"+map+"/Alignment/Final_Bam/"+id+"_SE_merged_"+map+".bam"
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
  
            for map in mappers:
                genom_cov={}#sample with genome coverage >40
                read_map={}#sample with mapped reads >90
                if not os.path.isdir(self.results+"/"+map+"/statistics/final_stat"):
                    os.makedirs(self.results+"/"+map+"/statistics/final_stat")   
                #check if there are merged Bam files or not and call the corresponding folder    
                if os.path.isdir(self.results+"/"+map+"/Alignment/Final_Bam"):
		    Final=self.results+"/"+map+"/Alignment/Final_Bam"
		else:
		    Final=self.results+"/"+map+"/Alignment/		            
		    for fl in os.listdir(Final):
			os.chdir(Final)	  
                    #if fl.endswith("_merged_bwa.bam") or fl.endswith("_final_bwa.bam")      
                        input=Final+"/"+fl
                        for ID in Sample:
                                        output1=self.results+"/"+map+"/statistics/"+ID+"_bwa_genomecov_bed.txt"
                                        #output2="+map+"_map_output+"/statistics/"+ID+"_bwa_genomecov>40.txt"
                                        output2=self.results+"/"+map+"/statistics/"+ID+"_bwa_genomecov=0.txt"
                                        os.system("""{}bedtools genomecov -ibam -bga {} -g {}.fasta > {}""".format(self.bedtools,input,self.ref_genome,output1)
                                        cov=os.system("awk '{ print $4 }' "+output1)
                                        genom_cov[ID]=cov          
                                        #os.system("""awk 'NF && $4>40' {} > {}""".format(output1,output2))
                                        os.system("""awk 'NF && $4<2' {} > {}""".format(output1,output3))
                                        output11=self."+map+"_map_output+"/statistics/"+ID+"_stats_bwa.txt"
                                        os.system("""{}samtools flagstat {} > {}""".format(self.samtools,input,output11))
                                        for fl in os.listdir(self.results+"/"+map+"/statistics/"): 
                                              if fl.endwith("_stats_"+map+".txt"):          
                                                  for line in open(fl).readlines():
                                                      line=line.strip()
                                                      if "mapped (" in line:
                                                              x=line.split("mapped")[1].split(":")[0].split("(")[1].split("%")[0]
                                                              y=float(x)
                                                              read_map[ID]=y                                                          
                                          for (i, j), (k, v) in product(genom_cov.items(), read_map.items()):   
                                                        if i==k and j>=90 or v>=40:
                                                            os.system("""cp {} {}""".format(self.results+"/"+map+"/statistics/"+ID+"_stats_"+map+".txt",self.results+"/"+map+"/statistics/final_stat"+ID+"_stats_"+map+".txt"))   
            return " "
    
    # Joint variant call SNP and indels using GATK
    def joint_variant_calling(self):
           print ("""
                  ==================================================================
	            Joint variant calling (SNP and indels) using GATK is starting !
	          ==================================================================
            """)				        
            for map in mappers:
                    if not os.path.isdir(self.results+"/"+map+"/Joint_Variants"):
                            os.makedirs(self.results+"/"+map+"/Joint_Variants")
		    if os.path.isdir(self.results+"/"+map+"/Alignment/Final_Bam"):
		    	Final=self.results+"/"+map+"/Alignment/Final_Bam"
		    else:
		        Final=self.results+"/"+map+"/Alignment/		            
		        for fl in os.listdir(Final):
			   os.chdir(Final)	       
                           input=Final+"/"+fl			  
                           for id in sample:
                                output=self.results+"/"+map+"/Joint_variants/"+id+"_GATK_snps_indels.vcf" 
                                os.system("""java -jar {}GenomeAnalysisTK.jar -R {}.fasta -T HaplotypeCaller -I {} --dbsnp {} -stand_call_conf 30 -stand_emit_conf 10.0 -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output))


            return "" 
    
    #structural variant identification(big deletion) using DELLY						  
    def structural_variation_calling(self):
	print ("""
                  ================================================
	             structural_variation_callingis starting !
	          ================================================
        """)				    
               
        for map in mappers:
                    if not os.path.isdir(self.results+"/"+map+"/struc_Variants"):
                            os.makedirs(self.results+"/"+map+"/struct_Variants")
                    if os.path.isdir(self.results+"/"+map+"/Alignment/Final_Bam"):
		    	Final=self.results+"/"+map+"/Alignment/Final_Bam"
		    else:
		        Final=self.results+"/"+map+"/Alignment/		            
		        for fl in os.listdir(Final):
			   os.chdir(Final)	       
                           input=Final+"/"+fl
                           for id in sample:
                                output=self.results+"/"+map+"/struc_Variants"+id+"_bwa_SV"
                                os.system("""{}delly call -t DEL -g {} -o {}.bcf {}""".format(self.DElly,self.ref_genome,output,input)
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
		for map in mappers:        
		    for fl in os.listdir(self.results+"/"+map+"/Joint_Variants/"):
			os.chdir(self.results+"/"+map+"/Joint_Variants/")
			file_name=os.path.splitext(fl)[0]
			input=fl
			output=self.results+"/"+map+"/Joint_Variants/"+file_name+".avinput"
			os.system("""{}convert2annovar.pl -format vcf4 {}} > {}.avinput""".format(self.annovar,input,output)
			os.system("""python -c "import sys;lines=sys.stdin.read();print lines.replace('gi|444893469|emb|AL123456.3|','Chromosome')" < {} > """ .format(output)           
		#run annovar gene annotation():
		    if not os.path.isdir(self.results+"/"+map+"/Annotation"):
			os.makedirs(self.results+"/"+map+"/Annotation")
		    for fl in os.listdir(self.results+"/"+map+"/Joint_Variants/"):
			if fl.endswith(".avinput"):
				os.chdir(self.results+"/"+map+"/Joint_Variants/")
				for id in sample:
					input=self.results+"/"+map+"/Joint_Variants/"+id+".vcf"
					output=self.results+"/"+map+"/Joint_Variants_annovar/"+id+".annovar"
					os.system("""perl {}table_annovar.pl --vcfinput {} {} --outfile {} --buildver H37Rv --nastring . --protocol refGene --operation g  --vcfinput""".format(self.annovar,input,self.TBdb,output)
		#annotate_variation.pl -filter -dbtype generic -genericdbfile h34rv_generic.txt -build hg19 -out ex1 h37rv.avinput example/        

		return""

if __name__ == '__main__':
    
    pipeline=BSMTBGan(global_Dir,user_ref,fastq_Dir,File_table,mappers)                                      
    print pipeline.spolygotyping()
    print pipeline.process_ref_genome()
    print pipeline.trimming()
    print pipeline.Mapping()
    print pipeline.realignment()
    print pipeline.base_qual_recal()
    print pipeline.PCR_dup_mark_remov()
    print pipeline.merge_bam()
    print pipeline.genome_coverag_mappability_stat()
    print pipeline.joint_variant_calling()
    print pipeline.structural_variation_calling()                                      
    print pipeline.Gene_annotation_annovar()

end_time = datetime.now()
new.write('Pipeline execution time: {}'.format(end_time - start_time))
new.close()
