#!/usr/bin/python
import os, glob
from shutil import copyfile

from pycompss.api.task import task
from pycompss.api.parameter import *
@task(fileName=FILE,binName=FILE, contigsFile=FILE_OUT, orfsFile=FILE_OUT, blastFile=FILE_OUT, kaijuFile=FILE_OUT)
def assembly_annotation(suffix,prefix,fileName,binName,contigsFile,orfsFile,blastFile,kaijuFile,cpus,assembler):
    
    #Make temp directory 
    directory= '~/DATMA_run_'+suffix
    cmd='mkdir '+directory
    print cmd
    os.system(cmd) 

    typeReads=prefix

    #Assembly process
    inputFile=binName
    outputFile=directory+'/contigs'
    assemblyOut=assembly(cpus,assembler,typeReads,inputFile,outputFile)
    
    cmd='cat '+outputFile+'/'+assemblyOut+' > '+contigsFile
    os.system(cmd) 

    #Open reading frames compute
    inputFile=outputFile+'/'+assemblyOut
    outputName=directory+'/orfs'
    orfs(inputFile,outputName)
    
    cmd='cat '+outputName+'.faa'+' > '+orfsFile
    os.system(cmd) 

    #nucleotide annotation
    inputFile=outputFile+'/'+assemblyOut
    outputName=directory+'/blastn'
    annotate_blastn(cpus,'/home/db/NT/NT_May_2017/nt',inputFile,outputName)
    
    cmd='cat '+outputName+' > '+blastFile
    os.system(cmd) 
    
    #Kaiju new method
    inputFile=outputFile+'/'+assemblyOut
    outputName=directory+'/kaiju'
    annotate_kaiju(cpus,'/home/users/andresb/tools/kaiju-master/kaijudb',inputFile,outputName)
    
    kronaName=outputName.split('.')
    kronaName=kronaName[0]+'.krona'
    
    cmd='cat '+kronaName+' > '+kaijuFile
    os.system(cmd) 
    
    
    #remove temp directory 
    directory= '~/DATMA_run_'+suffix
    cmd='rm -r '+directory
    os.system(cmd) 
    
# assembler
def assembly(cpus,assemblyTool,typeReads,inputName,outputName):
    if (assemblyTool == 'megahit'):	
        assemblyOut='final.contigs.fa'
        cmd='megahit -r '+inputName+' -o '+outputName+' -t '+cpus
        print cmd
        os.system(cmd)
    elif (assemblyTool == 'velvet'):
        assemblyOut='contigs.fa'
        cmd='velveth '+outputName+' 31 -'+typeReads+' -short '+inputName
        print cmd
        os.system(cmd)
        cmd='velvetg '+outputName+' -exp_cov auto -ins_length 260'
        print cmd
    	os.system(cmd)
    elif (assemblyTool == 'newbler'):
        assemblyOut='454LargeContigs.fna'
        cmd='runAssembly -mi 90 -ml 60 -cpu '+cpus+' -o '+outputName+' '+inputName
    	print cmd
    	os.system(cmd) 
    elif (assemblyTool == 'spades'):
        assemblyOut='contigs.fasta'
	cmd='cp '+inputName+' '+outputName+'_temp.'+typeReads
        print cmd
        os.system(cmd)

 	#cmd='spades.py --only-assembler -t '+cpus+' -o '+outputName+' -s '+outputName+'_temp.'+typeReads
        cmd='spades.py  -t '+cpus+' -o '+outputName+' -s '+outputName+'_temp.'+typeReads
        print cmd
        os.system(cmd)
    else:
        print "Error assembly no assigned"
        sys.exit(0)
        
    return assemblyOut

    print "PASS 4: assembling ... DONE"


#Open reading frames compute
def orfs(inputName,outputName):
    cmd='prodigal -i '+inputName+' -o '+outputName+' -a '+outputName+'.faa'
    print cmd
    os.system(cmd) 
    print "PASS 5: ORF ... DONE"

    
    
#Annotation for each bin using blastn
def annotate_blastn(cpus,database,inputName,outputName):
    cmd='blastn -query '+inputName+' -out '+outputName+' -db  '+database+' -outfmt 6 -num_alignments 3 -num_threads '+cpus
    print cmd
    os.system(cmd)
    print "PASS 7: Blastn ... DONE"

def annotate_kaiju(cpus,database,inputName,outputName):
    cmd='kaiju -z '+cpus+' -t '+database+'/nodes.dmp -f '+database+'/kaiju_db.fmi -i '+inputName+' -o '+outputName
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
    
    
    
    
    
