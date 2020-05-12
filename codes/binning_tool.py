#!/usr/bin/python
import os,glob
from map_genFm9 import *
from removal_tool import recoverySeq

#3.Binning the clean reads <CLAME>
def BinReads(direct,param,outName,flags,roundNum):
    cpus=str(param.cpus)
    typeReads=param.typeReads
    base=str(param.base)
    sizeBin=str(param.sizeBin)
    ld=str(param.ld)
    w=str(param.w)
    e=str(param.nu)
    clame_aux=str(param.clame_aux)
    
    directory=direct.globalOutput+'/'+direct.rounds
    inputFile=direct.globalOutput+'/'+outName.balance+'.'+typeReads
    outputFile=directory+'/'+outName.clame
    
    if flags.clame:
        if typeReads=='fastq' or typeReads=='sff':
            typeReads='fastq'
            formatReads='-fastq'
        else:
            typeReads='fasta'
            formatReads=' '
        
       
        #mapping 
        list2Exclude=direct.globalOutput+'/'+outName.balance+'.exclude'
        fm9ref=direct.globalOutput+'/'+outName.balance#+'.fm9'
        index=direct.globalOutput+'/'+outName.balance+'.index'
        numReads=sum(1 for line in open(index))
        totalFM=len(glob.glob(fm9ref+'*.fm9'))
        print("TOTAL FM9: "+str(totalFM)+" numReads: "+str(numReads))
        offset=param.block
        block=param.block
        mapping2FM9(cpus,base,numReads,inputFile,outputFile,typeReads,fm9ref,'-print',list2Exclude,clame_aux,totalFM,offset,block)
        
        #binning
        resultName=outputFile#+'.result'
        binning2FM9(cpus,resultName,outputFile,index,numReads,ld,sizeBin,clame_aux,list2Exclude)
        
        #excludeList =outputFile+'.exclude'
        #cmd='cat '+excludeList+' >>'+list2Exclude
        #print cmd
        #os.system(cmd) 
        
        
        #pathDataset = directory+'/'+outName.clame+'_*.list'  #input files
        #bins=glob.glob(pathDataset)
        #bins.sort(key=lambda x: int(x.split(outName.clame)[1].split(".")[0][1:]))
        #extractBinnedReads2(direct,param,outName,bins,flags,roundNum)

        #cmd='clame -b '+base+' -multiFasta '+inputFile+' -w '+w+' -ld '+ld+' -sizeBin '+sizeBin+' -e '+e+' -output '+outputFile+' -nt '+cpus+' '+formatReads+' '+clame_aux
        #print cmd
        #os.system(cmd) 
        
        #pathDataset = directory+'/'+outName.clame+'_*.fast*'  #input files
        #bins=glob.glob(pathDataset)
        #bins.sort(key=lambda x: int(x.split(outName.clame)[1].split(".")[0][1:]))
        #extractBinnedReads(direct,param,outName,bins,flags,roundNum)
        
        print("PASS 3: Binning ... DONE")
    else:
        outputFile+='_0.'+typeReads
        if not (os.path.exists(outputFile)):
            #symbolic link to the original file
            os.symlink(inputFile,outputFile)
        print('No binning stage selected')

def extractBinnedReads2(direct,param,outName,bins,flags,roundNum):
    cpus=param.cpus
    typeReads=param.typeReads
    
    directory= direct.globalOutput+'/'+direct.rounds
    inputFile=directory+'/'+outName.balance+'.'+typeReads
    outputFile=direct.globalOutput+'/'+outName.balance
    
    fsrc2=open(outputFile+'.list','w')     
    for binName in bins: #for each bin generated
        fsrc1 = open(binName,'r') 
        for line in fsrc1:
            fsrc2.write(line)
        fsrc1.close()
    fsrc2.close()
    
    recoverySeq(outputFile+'.list',inputFile,outputFile,typeReads,' ')
    
    #Remove temporal files
    cmd='wc -l '+outputFile+'.list'
    print(cmd)
    os.system(cmd) 
    
    cmd='rm '+outputFile+'.list'
    print(cmd)
    #os.system(cmd) 
    
def extractBinnedReads(direct,param,outName,bins,flags,roundNum):
    cpus=param.cpus
    typeReads=param.typeReads
    
    directory= direct.globalOutput+'/'+direct.rounds
    inputFile=directory+'/'+outName.balance+'.'+typeReads
    outputFile=direct.globalOutput+'/'+outName.balance
  
    firstC='>'
    if typeReads=='fasta':
        firstC='>'
        module=1
    elif typeReads=='fastq':
        firstC='@'
        module=4
    else:
        print('Quality Filter stage:ERROR IN THE CONFIGURATION FILE')
        print("Unknown reads' format")
        sys.exit(0)
        
    fsrc2=open(outputFile+'.list','w')     
    for binName in bins: #for each bin generated
        fsrc1 = open(binName,'r') 
        tipo=0
        for line in fsrc1:
            if line.startswith(firstC) and (tipo==0):
                fsrc2.write(line)
            tipo=(tipo+1) % module
        fsrc1.close()
    fsrc2.close()
    recoverySeq(outputFile+'.list',inputFile,outputFile,typeReads,' ')
    
    #Remove temporal files
    cmd='rm '+outputFile+'.list'
    print(cmd)
    os.system(cmd) 

