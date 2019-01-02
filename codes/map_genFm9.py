#!/usr/bin/python
import os,sys

#1.gen a FM9 file
def genFM9(cpus,inName,outName,typeReads):
    inputFile=inName
    outputFile=outName
    cmd=''
    #Generate the FM9 file
    if typeReads=='fastq' or typeReads=='sff':
        cmd='genFm9 -fastq -multiFasta '+inputFile+' -output '+outputFile+' -nt '+str(cpus)
    elif typeReads=='fasta':
        cmd='genFm9 -multiFasta '+inputFile+' -output '+outputFile+' -nt '+str(cpus)
    else:
        print "typeReads is not defined"
        sys.exit(0)

    print cmd
    os.system(cmd) 
    print "PASS : FM9 done"

def mapping2FM9(cpus,inName,outName,typeReads,fm9ref,print_enable):
    inputFile=inName
    outputFile=outName
    cmd=''
    #Mapping to FM9 file
    if typeReads=='fastq' or typeReads=='sff':
        cmd='mapping -fastq -b 20 -multiFasta '+inputFile+' -output '+outputFile+' -nt '+str(cpus)+' -fm9 '+fm9ref+' '+print_enable
    elif typeReads=='fasta':
        cmd='mapping -b 20 -multiFasta '+inputFile+' -output '+outputFile+' -nt '+str(cpus)+' -fm9 '+fm9ref+' '+print_enable
    else:
        print "typeReads is not defined"
        sys.exit(0)

    print cmd
    os.system(cmd) 
    print "PASS : mapping done"
