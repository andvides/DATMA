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
    forceMerge=False
    
    #16S-remove
    database_16s_fasta='~/DATMA/16sDatabases/ncbi/16SMicrobial.fasta'
    database_16s_fm9='~/DATMA/16sDatabases/ncbi/ncbi.fm9'
    RDP_path='~/DATMA/tools/RDPTools'
    aux_16s=" "
    useBWA=False
    
    #CLAME parameters
    bases='70,60,50,40,30,20'
    base=70
    ld=2
    sizeBin=2000
    nu=3 
    w=0
    clame_aux=" "
    fm9="unassigned"
    block=1
    
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
        elif line.startswith('-forceMerge'):
            param.forceMerge=True
            
        #16S parameters
        elif line.startswith('-useBWA'):
            param.useBWA=True
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
        elif line.startswith('-ld'):
            param.ld=int(words[1])
        elif line.startswith('-window'):
            param.w=int(words[1])
        elif line.startswith('-fm9'):
            param.fm9=words[1]
        elif line.startswith('-nu'):
            param.nu=float(words[1])
        elif line.startswith('-block'):
            param.block=int(words[1])
        elif line.startswith('-clame_aux'):
            param.clame_aux=" ".join(str(x) for x in words[1:])
            
        #BLAST parameters
        elif line.startswith('-database_nt'):
            param.database_nt=words[1]
        elif line.startswith('-database_nr'):
            param.database_nr=words[1]
        elif line.startswith('-blast_aux'):
            param.blast_aux=" ".join(str(x) for x in words[1:])
            
        #quast
        elif line.startswith('-use_ref'):
            param.use_ref=True
            param.ref_name=words[1]
        elif line.startswith('-quast_aux'):
            param.quast_aux=" ".join(str(x) for x in words[1:])
        elif line.startswith('-checkm_aux'): 
            param.checkm_aux=" ".join(str(x) for x in words[1:])
                
        #Assembler parameters
        elif line.startswith('-assembly'):
            param.assembly=words[1]
        elif line.startswith('-asm_aux'):
            param.asm_aux=words[1]
        elif line.startswith('-ORF_aux'):
            param.ORF_aux=" ".join(str(x) for x in words[1:])
        #Kaiju parameters
        elif line.startswith('-database_kaiju'):
            param.database_kaiju=words[1]
        elif line.startswith('-kaiju_aux'):
            param.kaiju_aux=" ".join(str(x) for x in words[1:])
            
    fsrc1.close()

    if error: #NO madatory inputs 
        print "ERROR IN THE CONFIGURATION FILE"
        sys.exit(0)
    else:
        print "PASS: Configuraion file loaded successfully"

        
        
#Update pipeline stages execution   
def pipelineFlow(param,stageFlags):
    startStage=param.startStage
    totalStages=param.totalStages
    
    if startStage<totalStages:
        if startStage==1:                   #full
            stageFlags.builddir=True
            stageFlags.quality=True
            stageFlags.seq_16sSeqs=True
            stageFlags.clame=True
            stageFlags.assemble=True            
            stageFlags.blastn=True            
            stageFlags.prodigal=True 
        elif startStage==2:                 #from 16s stage
            stageFlags.seq_16sSeqs=True
            stageFlags.clame=True
            stageFlags.assemble=True  
            stageFlags.blastn=True            
            stageFlags.prodigal=True 
        elif startStage==3:                 #from binning stage   
            stageFlags.clame=True
            stageFlags.assemble=True  
            stageFlags.blastn=True            
            stageFlags.prodigal=True 
        elif startStage==4:                 #from assemble    
            stageFlags.assemble=True  
            stageFlags.blastn=True            
            stageFlags.prodigal=True
        elif startStage==5:                 #from blast
            stageFlags.blastn=True            
            stageFlags.prodigal=True
        elif startStage==6:                 #from prodigal
            stageFlags.prodigal=True
        
        print 'STAGES: ',
        attrs = vars(stageFlags)
        print ', '.join("%s: %s" % item for item in attrs.items())

    else:
        print "ERROR IN THE CONFIGURATION FILE"
        sys.exit(0)   
        

#Buil directory
def builDirectory(direct,param,flags):
    if flags.builddir:
        startStage=param.startStage
        inputFile=param.inputFile
        direct.globalOutput=param.outputFile
        directory=direct.globalOutput
        if(os.path.exists(directory)):
            print 'Using an existing directory: ',directory 
            print 'start in: ',startStage,' Building: ',directory,' for: ',inputFile
            cmd="mkdir -p "+directory
            print cmd
            os.system(cmd)   
        else:
            print 'start in: ',startStage,' Building: ',directory,' for: ',inputFile
            cmd="mkdir "+directory
            print cmd
            os.system(cmd)
    else:
        direct.globalOutput=param.outputFile
        directory=direct.globalOutput
        print 'using: ',directory 

#Generate the final report        
def finalReport(binName):
    cmd="mkdir "+binName
    print cmd
    os.system(cmd)
    
    cmd="mv "+binName+'.* '+binName
    print cmd
    os.system(cmd)

def tablaReport(section):
    if (section=='header'):
        text="<!DOCTYPE html>\n<html>\n<head>\n<style>\ntable, th, td {\n\t border: 1px solid black;\n\tborder-collapse: collapse;\n}\n</style>\n</head>\n<body>\n"
    elif (section=='title'):
        text="<table style='width:70%'>\n<tr>\n\t<th>binNum</th>\n\t<th>binSize</th>\n\t<th>contigs</th>\n\t<th>genome</th>\n\t<th>genome</th>\n</tr>\n"
    elif(section=='end'):
        text="</table>\n</body>\n</html>\n"
    return text
    
        
def joinHtml():
    text='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"\n'
    text+='<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">\n'
    text+='<head>\n\t<meta charset="utf-8"/>\n</head>\n'
    text+='<body>\n'
    text+='<section>\n\t'
    text+='<h1>Bins report</h1>\n\t'
    text+='<iframe src="report.html" width="100%" frameborder="0" height="20%" ></iframe>\n'
    text+='</section>\n'
    text+='</section>\n'
    text+='<h1>Blastn report</h1>\n\t'
    text+='<iframe src="binsBlastn.html" width="100%" frameborder="0" height="600"></iframe>\n'
    text+='</section>\n'
    text+='</section>\n\t'
    text+='<h1>Kaiju report</h1>\n\t'
    text+='<iframe src="binsKaiju.html" width="100%" frameborder="0" height="600"></iframe>\n'
    text+="</section>\n</body>\n</html>\n"
    return text

def list2rawReads(binName,param,directory):
    typeReads=param.typeReads
    suffix=binName.split('.')
    ext=suffix[1]
    outFile=suffix[0]
    if typeReads=='illumina':
        input1=directory+"/clean/clean_reads_1.fastq"
        input2=directory+"/clean/clean_reads_1.fastq"
        cmd='selectFasta -fasta_sel -list '+binName+' -fastq '+input1+' >'+outFile+'_1.fastq'
        print cmd
        os.system(cmd)

        cmd='selectFasta -fasta_sel -list '+binName+' -fastq '+input2+' >'+outFile+'_2.fastq'
        print cmd
        os.system(cmd)
    elif typeReads=='fasta':
        input1=directory+"/clean/clean_reads.fastq"
        cmd='selectFasta -fasta_sel -list '+binName+' -fasta '+input1+' >'+outFile+'.fasta'
        print cmd
        os.system(cmd) 
        
    else:
        input1=directory+"/clean/clean_reads.fastq"
        cmd='selectFasta -fasta_sel -list '+binName+' -fastq '+input1+' >'+outFile+'.fastq'
        print cmd
        os.system(cmd) 


