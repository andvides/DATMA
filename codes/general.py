#!/usr/bin/python
import os,sys,subprocess

#Class to define the default configuration for the tool
#see the DATMA user manual for details
class parameters:
    #Basic configuration
    manual = 'save the basic inputs'
    totalStages=8
    startStage=1 
    endWorkflow=100
    inputFile= ''
    outputFile='output'
    cpus=1
    typeReads='sff'

    #Quality Control
    trimmomatic_path='~/DATMA/tools/Trimmomatic-0.38/'
    cleanTool='rapifilt'
    te=0    
    tb=0    
    lq=30   
    rq=30   
    m=70    
    wq=2    
    quality_aux=" "
    
    #Flash
    fb=5
    flash_aux=" "
    
    #16S-remove
    database_16s_fasta='~/DATMA/16sDatabases/ncbi/16SMicrobial.fasta'
    database_16s_fm9='~/DATMA/16sDatabases/ncbi/ncbi.fm9'
    RDP_path='~/DATMA/tools/RDPTools'
    aux_16s=" "
    
    #CLAME parameters
    bases='70,60,50,40,30,20'
    base=70
    ld=2
    sizeBin=2000
    nu=3 
    w=0
    clame_aux=" "
    
    #Assembly options
    assembly='spades'	    
    assemblyOut='contigs.fa'    
    asm_aux=" "
    ORF_aux=" "
    
    #quast
    use_ref=False
    ref_name=' '
    quast_aux=" "
    checkm_aux=" "

    #BLAST
    database_nt='~/DATMA/blastdb/nt'
    database_nr='~/DATMA/blastdb/nr'
    blast_aux=" "
    
    #Kaiju
    database_kaiju='~/DATMA/tools/kaiju/kaijudb'
    kaiju_aux=" "

#Tools that compose the DATMA framework  
#all the tools need to be added to user PATH
class tools:
    manual = 'save the tools names'
    names=['selectFasta','rapifilt','flash',
           'bwa','samtools','mapping','genFm9','clame',
           'megahit','spades.py','velvetg',
           'blastn','kaiju','quast.py','prodigal',
           'ktImportBLAST','checkm']

#Directorie path for the outputs
class directories:
    globalOutput='output'
    filterQuality='clean'
    removal='16sSeq'
    rounds='round_'
    blast='blast'
    orfs='prodigal'

#Output names for the files generates    
class outputNames:
    manual = 'save the out names for each output tool'
    clean='clean_reads'   
    seq_16sSeqs='all_16S'
    clame= 'clameBin'
    assemble='contigs'
    blastn='taxo.tab' 
    prodigal='contigs_orfs'
    balance='readsForbin'
    kaiju='taxo.kaiju'

#flags to enable or disable the DATMA stages    
class stageFlags:
    manual = 'save the out execution flag for each output tool'
    builddir=False
    quality=False
    seq_16sSeqs=False
    clame=False
    assemble=False
    blastn=False 
    prodigal=False


# check install tools
def is_tool(tool):
    for items in tool.names:
        rc = subprocess.call(['which', items])
        if rc:
            print items,' missing in path!'
            sys.exit(0)
    print 'PASS: all the tools are installed'      
    
# Read configuration file
def readConfigFile(fileName, param):
    error=True
    fsrc1 = open(fileName,'r') 
    for line in fsrc1:
        words=line.split("\n") 
        words=words[0] 
        words=words.split()
        #Basic parameters
        if line.startswith('-start_in'):
            initial_stage=int(words[1])
            param.startStage=initial_stage
        elif line.startswith('-end_in'):
            param.endWorkflow=int(words[1])
        elif line.startswith('-inputFile'): #mandatory!!!
            error=False
            inputFile=words[1]
            param.inputFile=inputFile
        elif line.startswith('-outputDir'):
            directory=words[1]
            param.outputFile=directory
        elif line.startswith('-cpus'):
            cpus=words[1]
            param.cpus=cpus
        elif line.startswith('-typeReads'): 
            typeReads=words[1]
            param.typeReads=typeReads 
            if not (param.typeReads=='sff' or param.typeReads=='illumina' or param.typeReads=='fastq' or param.typeReads=='fasta'):
                print "ERROR IN THE CONFIGURATION FILE"
                print "unsupported reads' format"
                sys.exit(0)
        
        #Quality control tools
        elif line.startswith('-trimmomatic_path'):
            param.trimmomatic_path=words[1]
        elif line.startswith('-cleanTool'):
            param.cleanTool=words[1]    
        elif line.startswith('-te'):
            param.te=int(words[1])
        elif line.startswith('-tb'):
            param.tb=int(words[1])
        elif line.startswith('-lq'):
            param.lq=int(words[1])
        elif line.startswith('-rq'):
            param.rq=int(words[1])
        elif line.startswith('-m'):
            param.m=int(words[1])
        elif line.startswith('-wq'):
            param.wq=int(words[1])
        elif line.startswith('-quality_aux'):
            param.quality_aux=" ".join(str(x) for x in words[1:])

        #Flash
        elif line.startswith('-fb'):
            param.fb=int(words[1]) 
        elif line.startswith('-flash_aux'):
            param.flash_aux=" ".join(str(x) for x in words[1:])           
        
        #16S parameters
        elif line.startswith('-database_16s_fasta'):
            param.database_16s_fasta=words[1]
        elif line.startswith('-database_16s_fm9'):
            param.database_16s_fm9=words[1]
        elif line.startswith('-RDP_path'):
            param.RDP_path=words[1]
        elif line.startswith('-aux_16s'):
            param.aux_16s=" ".join(str(x) for x in words[1:])
            
        #CLAME parameters        
        elif line.startswith('-bases'):            
            param.bases=str(words[1])
        elif line.startswith('-sizeBin'):            
            sizeBin=int(words[1])
            param.sizeBin=sizeBin
        elif line.startswi
