#!/usr/bin/python
import os,sys

def merge(param,inName,outName,directory):
    #merge files
    cpus=str(param.cpus)
    mf=str(param.fb)
    aux_flash=str(param.flash_aux);
    cmd="flash -q -m "+mf+" "+inName+"_1.fastq "+inName+"_2.fastq  -t "+cpus+' -o '+outName+' -d '+directory+' '+aux_flash
    print cmd
    os.system(cmd) 
    print "***********************************"
    print "PASS: Merge Files  ... DONE"
    print "***********************************"
    
