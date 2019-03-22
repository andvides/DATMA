#!/usr/bin/python
import os, glob
from shutil import copyfile

from pycompss.api.task import task
from pycompss.api.parameter import *
@task(binName=FILE, contigsFile=FILE_OUT, contigsFile_stat=FILE_OUT, orfsFile=FILE_OUT, blastFile=FILE_OUT, kaijuFile=FILE_OUT)
def assembly_annotation(suffix,prefix,binName,contigsFile,contigsFile_stat,orfsFile,blastFile,database_nt,kaijuFile,database_kaiju,cpus,assembler,use_ref,ref_name,asm_aux,ORF_aux,blast_aux,kaiju_aux,quast_aux):
    
    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    #Edit these lines if the database_nt and database_kaiju directories are different in the workerers
    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    database_nt_worker=database_nt;#full path to the nt database
    database_kaiju_worker=database_kaiju;#full path to the kaiju database
    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    
    #Make temp directory 
    directory= '~/DATMA_run_'+suffix
    cmd='mkdir '+directory
    print cmd
    os.system(cmd) 

    typeReads=prefix

    #Assembly process
    inputFile=binName
    outputFile=directory+'/contigs'
    assemblyOut=assembly(cpus,assembler,typeReads,inputFile,outputFile,asm_aux)
    
    cmd='cat '+outputFile+'/'+assemblyOut+' > '+contigsFile
    os.system(cmd) 

    #assembler quality
    inputFile=outputFile+'/'+assemblyOut
    dir_out=directory+'/quast'
    ass_quast(inputFile,dir_out,use_ref,ref_name,quast_aux)
    
    cmd='cat '+dir_out+'/report.html > '+contigsFile_stat
    os.system(cmd) 

    #Open reading frames compute
    inputFile=outputFile+'/'+assemblyOut
    outputName=directory+'/orfs'
    orfs(inputFile,outputName,ORF_aux)
    
    cmd='cat '+outputName+'.faa'+' > '+orfsFile
    os.system(cmd) 

    #nucleotide annotation
    inputFile=outputFile+'/'+assemblyOut
    outputName=directory+'/blastn'
    annotate_blastn(cpus,database_nt_worker,inputFile,outputName,blast_aux)
    
    cmd='cat '+outputName+' > '+blastFile
    os.system(cmd) 
    
    #Kaiju new method
    inputFile=outputFile+'/'+assemblyOut
    outputName=directory+'/kaiju'
    annotate_kaiju(cpus,database_kaiju_worker,inputFile,outputName,kaiju_aux)
    
    kronaName=outputName.split('.')
    kronaName=kronaName[0]+'.krona'
    
    cmd='cat '+kronaName+' > '+kaijuFile
    os.system(cmd) 
    
    
    #remove temp directory 
    directory= '~/DATMA_run_'+suffix
    cmd='rm -r '+directory
    os.system(cmd) 
    
# assembler
def assembly(cpus,assemblyTool,typeReads,inputName,outputName,asm_aux):
    if (assemblyTool == 'megahit'):	
        assemblyOut='final.contigs.fa'
        cmd='megahit -r '+inputName+' -o '+outputName+' -t '+cpus+' '+asm_aux
        print cmd
        os.system(cmd)
    elif (assemblyTool == 'velvet'):
        assemblyOut='contigs.fa'
        cmd='velveth '+outputName+' 31 -'+typeReads+' -short '+inputName+' '+asm_aux
        print cmd
        os.system(cmd)
        cmd='velvetg '+outputName+' -exp_cov auto -ins_length 260'+' '+asm_aux
        print cmd
    	os.system(cmd)
    elif (assemblyTool == 'newbler'):
        assemblyOut='454LargeContigs.fna'
        cmd='runAssembly -mi 90 -ml 60 -cpu '+cpus+' -o '+outputName+' '+inputName+' '+asm_aux
    	print cmd
    	os.system(cmd) 
    elif (assemblyTool == 'spades'):
        assemblyOut='contigs.fasta'
        cmd='cp '+inputName+' '+outputName+'_temp.'+typeReads
        print cmd
        os.system(cmd)

        #cmd='spades.py --only-assembler -t '+cpus+' -o '+outputName+' -s '+outputName+'_temp.'+typeReads
        cmd='spades.py  -t '+cpus+' -o '+outputName+' -s '+outputName+'_temp.'+typeReads+' '+asm_aux
        print cmd
        os.system(cmd)
    else:
        print "Error assembly no assigned"
        sys.exit(0)
        
    return assemblyOut

    print "PASS 4: assembling ... DONE"

#assembler quality
def ass_quast(cotigs_in,dir_out,ref,ref_name,quast_aux):
    if ref:
        cmd='quast.py -o '+dir_out+' '+cotigs_in+' -r '+ref_name+' '+quast_aux 
        print cmd
        os.system(cmd)
    else:
        cmd='quast.py -o '+dir_out+' '+cotigs_in+' '+quast_aux
        print cmd
        os.system(cmd)

#Open reading frames compute
def orfs(inputName,outputName,ORF_aux):
    cmd='prodigal -i '+inputName+' -o '+outputName+' -a '+outputName+'.faa'+' '+ORF_aux
    print cmd
    os.system(cmd) 
    print "PASS 5: ORF ... DONE"

    
    
#Annotation for each bin using blastn
def annotate_blastn(cpus,database,inputName,outputName,blast_aux):
    cmd='blastn -query '+inputName+' -out '+outputName+' -db  '+database+' -outfmt 6 -num_alignments 3 -num_threads '+cpus+' '+blast_aux
    print cmd
    os.system(cmd)
    print "PASS 7: Blastn ... DONE"

def annotate_kaiju(cpus,database,inputName,outputName,kaiju_aux):
    cmd='kaiju -z '+cpus+' -t '+database+'/nodes.dmp -f '+database+'/kaiju_db.fmi -i '+inputName+' -o '+outputName+' '+kaiju_aux
    print cmd
    os.system(cmd)
        
    kronaName=outputName.split('.')
    kronaName=kronaName[0]+'.krona'
    cmd='kaiju2krona -t '+database+'/nodes.dmp -n '+database+'/names.dmp -i '+outputName+' -o '+kronaName
    print cmd
    os.system(cmd)    

def makeReport(typeReads,inputName):
    char='>'
    if (typeReads=='fastq'):
        char='@'
    count=0
    bases=0
    enable=True
    
    if (os.path.exists(inputName)):
        fsrc1 = open(inputName,'r') 
        for line in fsrc1:
            if enable:
                if(line.startswith(char)): #fasta or fastq
                    count+=1
                elif(line.startswith('+')): #fastq
                    enable=False
                else:
                    words=line.split("\n") 
                    words=words[0] 
                    bases+=len(words)
            else: #quality
                enable=True
        fsrc1.close()
        return str(count),str(bases)
    else:
        return str('Error'),str(inputName)
    
