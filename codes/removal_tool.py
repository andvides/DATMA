#!/usr/bin/python
import os,sys
from map_genFm9 import *

#STAGE2.16Sremoval from the clean reads 
def seq_16Sremoval(direct,param,outName,flags):
    database_fasta=param.database_16s_fasta
    database_fm9=param.database_16s_fm9
    aux_16s=str(param.aux_16s)
    
    directory= direct.globalOutput
    inputFile=directory+'/'+direct.filterQuality+'/'+outName.clean
    outFile=directory+'/'+direct.removal+'/'+outName.seq_16sSeqs
   
    cpus=param.cpus
    typeReads=param.typeReads
    assemblyTool=param.assembly
    RDP_path=param.RDP_path
    useBWA=param.useBWA
    clame_aux=param.clame_aux

    if flags.seq_16sSeqs:
        #make directory to save the results
        cmd="mkdir "+directory+'/'+direct.removal
        print(cmd)
        os.system(cmd)
        
        firstC='>'
        if typeReads=='fasta':
            firstC='>'
        elif typeReads=='fastq':
            firstC='@'
        else:
            print('Quality Filter stage:ERROR IN THE CONFIGURATION FILE')
            print("Unknown reads' format")
            sys.exit(0)
            
        #Mapping the clean reads against the 16S_database
        if useBWA:
            temp_BWA=directory+'/'+direct.removal+'/temp_BWA'
            cmd='bwa mem -t '+cpus+' '+database_fasta+' '+inputFile+"."+typeReads+' > '+temp_BWA+'.sam'+' '+aux_16s
            print(cmd)
            os.system(cmd)
        
            cmd='samtools view -S -F4 '+' '+temp_BWA+'.sam > '+temp_BWA+'.map'
            print(cmd)
            os.system(cmd)

            cmd="awk '{print $1}' "+temp_BWA+".map > "+outFile+'.list'
            print(cmd)
            os.system(cmd)

            cmd='rm '+temp_BWA+'*'
            print(cmd)
            os.system(cmd)
        else:
            list2Exclude=""
            mapping2FM9(cpus,20,inputFile+"."+typeReads,outFile+"_1",typeReads,database_fm9,'-print',list2Exclude,clame_aux)
            cmd= "awk '{print $1}' "+outFile+"_1.links | awk -F '"+firstC+"' '{print $2}' > "+outFile+'_1.list'
            print(cmd)
            os.system(cmd)   
        
            #Generate the FM9 from the clean reads
            genFM9(cpus,inputFile+"."+typeReads,inputFile,typeReads)
            cleanFM9=inputFile+".fm9"
            mapping2FM9(cpus,20,database_fasta,outFile+"_2",'fasta',cleanFM9,'-print',list2Exclude,clame_aux)
            decoseq(firstC,inputFile+'.index',outFile+"_2.result",outFile+"_2.list")

            #concatenate 16s sequences found
            cmd="cat "+outFile+"_1.list "+outFile+"_2.list | sed '/^\s*$/d' > "+outFile+".list"
            print(cmd)
            os.system(cmd)
            #cmd="cat "+outFile+"_1.list | sed '/^\s*$/d' > "+outFile+".list"

            #remove temp files
            cmd="rm "+outFile+"_*.*"
            print(cmd)
            os.system(cmd) 
        
            cmd="rm "+cleanFM9
            print(cmd)
            os.system(cmd) 
        
        
        #Recovery the sequences names
        recoverySeq(outFile+".list",inputFile+"."+typeReads,outFile,typeReads,"-fasta_sel")    

        #removingSeqs from the pool of the reads
        outputName=directory+'/'+outName.balance
        recoverySeq(outFile+".list",inputFile+"."+typeReads,outputName,typeReads," ")    

        #count number of 16S-sequences found
        if typeReads=='fastq':
            firstC='+'
        num_lines=0
        fsrc1 = open(outFile+'.'+typeReads,'r') 
        for line in fsrc1:
            if line.startswith(firstC):
                num_lines+=1
        fsrc1.close()
        
        #RDP classifier
        #RDP_path='/home/software/src/rdp_classifier_2.7/dist'
        inputName=outFile+'.'+typeReads
        outputName=outFile+'.rdp'
        if RDP_path=='none':
            if typeReads=='fastq': #new rdp supports only fasta files
                temp_name=outFile
                cmd='rapifilt -f -fastq '+inputName+' -o '+temp_name
                print("pass to fasta files")
                os.system(cmd) 
                cmd='rdp_classifier'+' -o '+outputName+' -q '+temp_name+'.fasta'
            else:
                cmd='rdp_classifier'+' -o '+outputName+' -q '+inputName
        else:
            cmd='java -Xmx1g -jar '+RDP_path+'/classifier.jar'+' -c 0.5 '+' -o '+outputName+' '+inputName
        print(cmd)
        os.system(cmd) 
        
        #html report
        #cmd='ktImportRDP '+outputName+' -o '+'16S_rdp.html'
        #print(cmd)
        #os.system(cmd) 
    
        print("PASS 2: "+str(num_lines)+" sequences found")
    else:
        inputFile+='.'+typeReads
        outputFile=directory+'/'+outName.balance+'.'+typeReads
        if not (os.path.exists(outputFile)):
            #symbolic link to the original file
            os.symlink(inputFile,outputFile)
        print('No 16S remove stage selected')

def recoverySeq(inputList,inputFile,outFile,typeReads,fasta_sel):        
    cmd="selectFasta "+fasta_sel+" -list "+inputList+" -"+typeReads+" "+inputFile+' >'+outFile+'.'+typeReads
    print(cmd)
    os.system(cmd) 

    
def decoseq(firstC,index,result,outFile):
    fsrc1 = open(result,'r') 
    results=set()
    for line in fsrc1:
        words=line.split('\n')
        words=words[0].split()
        for items in words[1:]:
            results.add(int(items))
    fsrc1.close()

    fsrc2 = open(outFile,'w')
    fsrc3 = open(index,'r') 
    for line in fsrc3:
        words=line.split('\n')
        words=words[0].split()
        name=words[0].split(firstC)
        name=name[1]
        index=int(words[1])
        if index in results:
            fsrc2.write(name+'\n')
    fsrc2.close()
    fsrc3.close()

