#!/usr/bin/python
import os,sys,glob

#1.gen a FM9 file
def genFM9(cpus,inName,outName,typeReads, block):
    inputFile=inName
    outputFile=outName
    cmd=''
    perblock=' ';
    if block>1:
        perblock=' -size '+str(block)
    #Generate the FM9 file
    if typeReads=='fastq' or typeReads=='sff':
        cmd='genFm9 -fastq -multiFasta '+inputFile+' -output '+outputFile+perblock
    elif typeReads=='fasta':
        cmd='genFm9 -multiFasta '+inputFile+' -output '+outputFile+perblock
    else:
        print("typeReads is not defined")
        sys.exit(0)

    print(cmd)
    os.system(cmd) 
    print("PASS : FM9 done")
    
def mapping2FM9(cpus,bases,numReads,inName,outName,typeReads,fm9ref,print_enable,list2Exclude,clame_aux,totalFM,offset,block):
    inputFile=inName
    outputFile=outName
    totalReads=str(numReads+1)
    cmd1=''
    exclude=' -list2Exclude '+list2Exclude 

    cmd1='mapping -b '+str(bases)+' -multiFasta '+inputFile+' -nt '+str(cpus)+' '+print_enable+exclude
    
    if typeReads=='fastq' or typeReads=='sff':
        cmd1+=' -fastq '
    else:
        print("typeReads is not defined")
        sys.exit(0)

    #Mapping to FM9 file
    if block>1:
        cmd1+=' -size '+str(block)+' -totalReads '+totalReads
        print("-----------------------------------------------------------------------")
        for i in range(totalFM):
            groupnum=int(i*offset)
            cmd=cmd1+' -offsetFM9 '+str(groupnum)+' -output '+outputFile+'_'+str(i)+' -fm9 '+(fm9ref+str(i)+'.fm9')
            print(cmd)
            print("-----------------------------------------------------------------------")
            os.system(cmd) 
    else:
        print("-----------------------------------------------------------------------")
        cmd=cmd1+' -output '+outputFile+' -fm9 '+(fm9ref+'.fm9') +' -totalReads '+totalReads
        print(cmd)
        print("-----------------------------------------------------------------------")
        os.system(cmd) 

    print("PASS : mapping done")

def binning2FM9(cpus,resultName,outName,index,numReads,ld,sizeBin,clame_aux,list2Exclude):
    totalReads=str(numReads+1)
    results=''
    for files in glob.glob(resultName+'*.resultb'):
        results+=files+','
    
    if not (results==''):
        cmd='binning -i '+index+' -n '+totalReads+' -nt '+cpus+' -ld '+str(ld)+' -rt '+results[:-1]+' -sizeBin '+str(sizeBin)+' -output '+outName+' '+clame_aux+' > '+outName+'.binning'
        print(cmd)
        os.system(cmd) 
        
        cmd='rm '+resultName+'*'
        #os.system(cmd) 
        print("PASS : binning done")
    
        excludeList =outName+'.exclude'
        cmd='cat '+excludeList+' >>'+list2Exclude
        print(cmd)
        os.system(cmd) 

