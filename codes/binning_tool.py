#!/usr/bin/python
import os,glob

#3.Binning the clean reads <CLAME>
def BinReads(direct,param,outName,flags,roundNum):
    cpus=str(param.cpus)
    typeReads=param.typeReads
    base=str(param.base)
    sizeBin=str(param.sizeBin)
    ld=str(param.ld)
    w=str(param.w)
    e=str(param.nu)
    
    directory=direct.globalOutput+'/'+direct.rounds
    inputFile=directory+'/'+outName.balance+'.'+typeReads
    outputFile=directory+'/'+outName.clame
    
    if flags.clame:
        if typeReads=='fastq' or typeReads=='sff':
            typeReads='fastq'
            formatReads='-fastq'
        else:
            typeReads='fasta'
            formatReads=' '
            
        cmd='clame -b '+base+' -multiFasta '+inputFile+' -w '+w+' -ld '+ld+' -sizeBin '+sizeBin+' -e '+e+' -output '+outputFile+' -nt '+cpus+' '+formatReads
        print cmd
        os.system(cmd) 
        
        pathDataset = directory+'/'+outName.clame+'_*.fast*'  #input files
        bins=glob.glob(pathDataset)
        bins.sort(key=lambda x: int(x.split(outName.clame)[1].split(".")[0][1:]))
        extractBinnedReads(direct,param,outName,bins,flags,roundNum)
        
        print "PASS 3: Binning ... DONE"
    else:
        outputFile+='_0.'+typeReads
        if not (os.path.exists(outputFile)):
            #symbolic link to the original file
            os.symlink(inputFile,outputFile)
        print 'No binning stage selected'

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
        print 'Quality Filter stage:ERROR IN THE CONFIGURATION FILE'
        print "Unknown reads' format"
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
    from removal_tool import recoverySeq
    recoverySeq(outputFile+'.list',inputFile,outputFile,typeReads,' ')
    
    #Remove temporal files
    cmd='rm '+outputFile+'.list'
    print cmd
    os.system(cmd) 
