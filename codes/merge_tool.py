#!/usr/bin/python
import os,sys

def merge(param,inName,outName,directory):
    #merge files
    cpus=str(param.cpus)
    mf=str(param.fb)
    cmd="flash2 -q -m "+mf+" "+inName+"_1.fastq "+inName+"_2.fastq  -t "+cpus+' -o '+outName+' -d '+directory
    print cmd
    os.system(cmd) 
    print "***********************************"
    print "PASS: Merge Files  ... DONE"
    print "***********************************"
    
