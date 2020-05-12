#!/usr/bin/python
import os,sys,argparse

class parameters:
    manual = 'save the basic inputs'
    endWorkflow=100
    totalStages=8
    inputFile= ''
    outputFile='output'
    cpus=1
    typeReads='sff'      
    combine=' '
    checkm_aux=''
    lineage=''
    delete=False;

def readConfigFile(fileName, param):
    error=True
    fsrc1 = open(fileName,'r') #abre el archivo 
    for line in fsrc1:
        words=line.split("\n") 
        words=words[0] 
        words=words.split()
        #basic parameters
        if line.startswith('-inputFile'): #mandatory!!!
            error=False
            inputFile=words[1]
            param.inputFile=inputFile
        elif line.startswith('-end_in'):
            param.endWorkflow=int(words[1])
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
                print("ERROR IN THE CONFIGURATION FILE")
                print("unsupported reads' format")
        elif line.startswith('-combine'): 
            param.combine='-c' 
        elif line.startswith('-checkm_aux'): 
            param.checkm_aux=words[1]            
        elif line.startswith('-lineage'): 
            param.lineage=words[1] 
        elif line.startswith('-delete'): 
            param.delete=True

def removeFiles(directory,typeReads):
    if (os.path.exists(directory)):
        #Clean
        #delete clean directory, only conserved clean_reads compressed
        #pathDataset=directory+'/clean/'
        #cmd='gzip '+pathDataset+'/clean_reads.'+typeReads
        #print(cmd)
        #os.system(cmd)    
                    
        pathDataset=directory+'/clean/'
        cmd='rm '+pathDataset+'/*.txt '+pathDataset+'/*.'+typeReads
        print(cmd)
        os.system(cmd)                  
        #delete 16s rRNA
        pathDataset=directory+'/16sSeq/'
        cmd='rm '+pathDataset+'/*.'+typeReads
        print(cmd)
        os.system(cmd) 
                    
        #readsForbin
        pathDataset=directory
        cmd='rm '+pathDataset+'/*.fm9 '+pathDataset+'/*.'+typeReads #pathDataset/*.index
        os.system(cmd)
                    
        #round
        pathDataset=directory+'/round_*'
        cmd='rm '+pathDataset+'/*.index '+pathDataset+'/*.baseIndex '+pathDataset+'/*.links '+pathDataset+'/*.resultb '+pathDataset+'/*.'+typeReads
        print(cmd)
        os.system(cmd)
                    
        #busco
        pathDataset=directory+'/busco'
        cmd='rm -rf '+pathDataset+'/run_*/'
        print(cmd)
        os.system(cmd)
                    
        #CheckM
        pathDataset=directory+'/checkm_out'
        cmd='rm -rf '+pathDataset+'/bins/ '+pathDataset+'/storage '+pathDataset+'/lineage.ms'
        print(cmd)
        os.system(cmd)
        
        
        
    else:
        print(directory+" directory does not exists")
    

        
def yes_or_no(question):
    question=question
    check = str(input(question+" (Y/N): ")).lower().strip()
    try:
        if check[0] == 'y':
            return True
        elif check[0] == 'n':
            return False
        else:
            print('Invalid Input')
            return yes_or_no(question)
    except Exception as error:
        print("Please enter valid inputs")
        print(error)
        return yes_or_no(question)
    
def yes_or_no2(question):
    question=question
    check = str(raw_input(question+" (Y/N): ")).lower().strip()
    try:
        if check[0] == 'y':
            return True
        elif check[0] == 'n':
            return False
        else:
            print('Invalid Input')
            return yes_or_no(question)
    except Exception as error:
        print("Please enter valid inputs")
        print(error)
        return yes_or_no(question)

if __name__ == "__main__":
    #Final report
    description='DATMA: Distributed AuTomatic Metagenomc Assembly and Annotation framework delete no essencial files script \n'
    epilog='Authors: '
    epilog+='Benavides A, Alzate JF and Cabarcas F \n'
    
    question="You will delete most intermediate files to save disk space.\n"
    question+="However, it would prevent DATMA from restarting from an intermediate point\n"
    question+="Are you sure you want to continue?"

    
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-f', '--file', help='configuration file', required=True)

    args = parser.parse_args()
    param=parameters()
    
    if args.file :
        readConfigFile(args.file,param)
        typeReads=param.typeReads
        if not (param.typeReads=='fastq' or param.typeReads=='fasta'):
            typeReads='fast*'
            
        directory=param.outputFile
        
        if sys.version_info.major == 3:
            if(yes_or_no(question)):
                removeFiles(directory,typeReads)
            else:
                print ("NO changes done")
        else:
            if(yes_or_no2(question)):
                removeFiles(directory,typeReads)
            else:
                print ("NO changes done")
        
        
    else:
        print("Usage deleteDATMAfiles -f configurationFile")
        



