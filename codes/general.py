#!/usr/bin/python
import os,sys,subprocess

#Class to define the default configuration for the tool
#see the DATMA user manual for details
class parameters:
    #Basic configuration
    manual = 'save the basic inputs'
    totalStages=8
    startStage=1 
    inputFile= ''
    outputFile='output'
    cpus=1
    typeReads='sff'

    #Quality Control
    cleanTool='rapifilt'
    te=0    
    tb=0    
    lq=30   
    rq=30   
    m=70    
    wq=2    

    #Flash
    fb=5
    
    #16S-remove
    database_16s_fasta='~/16sDatabase/rfam/RFAM_db.fasta'
    database_16s_fm9='~/16sDatabase/rfam/rfam.fm9'
    
    #CLAME parameters
    bases='70,60,50,40,30,20'
    base=70
    ld=2
    sizeBin=10000
    nu=3 
    w=0
    
    #Assembly options
    assembly='velvet'	    
    assemblyOut='contigs.fa'    
    
    #BLAST
    database_nt='~/nt'
    database_nr='~/nr'
    
    #Kaiju
    database_kaiju='~/kaijudb'
    

#Tools that compose the DATMA framework  
#all the tools need to be added to user PATH
class tools:
    manual = 'save the tools names'
    names=['selectFasta','rapifilt','mapping','genFm9','flash','clame','megahit','spades.py','velvetg','blastn','kaiju','prodigal','ktImportBLAST']

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
	
	#Flash
        elif line.startswith('-fb'):
            param.fb=int(words[1])                
        
        #16S parameters
        elif line.startswith('-database_16s_fasta'):
            param.database_16s_fasta=words[1]
        elif line.startswith('-database_16s_fm9'):
            param.database_16s_fm9=words[1]
        
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
        elif line.startswith('-nu'):
            param.nu=float(words[1])
            
        #BLAST parameters
        elif line.startswith('-database_nt'):
            param.database_nt=words[1]
        elif line.startswith('-database_nr'):
            param.database_nr=words[1]
        
        #Assembler parameters
        elif line.startswith('-assembly'):
            param.assembly=words[1]
        
        #Kaiju parameters
        elif line.startswith('-kd'):
            param.database_kaiju=words[1]
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
        if startStage==1:
            stageFlags.builddir=True
            stageFlags.quality=True
            stageFlags.seq_16sSeqs=True
            stageFlags.clame=True
            stageFlags.assemble=True            
            stageFlags.blastn=True            
            stageFlags.prodigal=True 
        elif startStage==2:
            stageFlags.seq_16sSeqs=True
            stageFlags.clame=True
            stageFlags.assemble=True  
            stageFlags.blastn=True            
            stageFlags.prodigal=True 
        elif startStage==3:   
            stageFlags.clame=True
            stageFlags.assemble=True  
            stageFlags.blastn=True            
            stageFlags.prodigal=True 
        elif startStage==4:    
            stageFlags.assemble=True  
            stageFlags.blastn=True            
            stageFlags.prodigal=True
        elif startStage==5:    
            stageFlags.blastn=True            
            stageFlags.prodigal=True
        elif startStage==6:
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
            print "ERROR IN THE CONFIGURATION FILE"
            print 'Using an existing directory: ',directory 
            sys.exit(0)   
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
