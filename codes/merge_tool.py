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
    
    if(param.forceMerge):    #Force not combined to merge
        cmd="mergeNotCombined "+directory+"/"+outName+".notCombined_1.fastq "+directory+"/"+outName+".notCombined_2.fastq "+"NNN "+"> "+directory+"/"+outName+".force.fastq"
        print cmd
        os.system(cmd) 
    
        #Mix the ouptuts
        cmd="cat "+directory+"/"+outName+".force.fastq >> "+directory+"/"+outName+'.extendedFrags.fastq'
        print cmd
        os.system(cmd) 
    
    print "***********************************"
    print "PASS: Merge Files  ... DONE"
    print "***********************************"
    
def mergeGeneric(param,inName1,inName2,outName,directory):
    #merge files
    cpus=str(param.cpus)
    mf=str(param.fb)
    aux_flash=str(param.flash_aux);
    cmd="flash -q -m "+mf+" "+inName1+" "+inName2+"  -t "+cpus+' -o '+outName+' -d '+directory+' '+aux_flash
    print cmd
    os.system(cmd) 
    
    if(param.forceMerge): #Force not combined to merge
        cmd="mergeNotCombined "+directory+"/"+outName+".notCombined_1.fastq "+directory+"/"+outName+".notCombined_2.fastq "+"NNN "+"> "+directory+"/"+outName+".force.fastq"
        print cmd
        os.system(cmd) 
    
        #Mix the ouptuts
        cmd="cat "+directory+"/"+outName+".force.fastq >> "+directory+"/"+outName+'.extendedFrags.fastq'
        print cmd
        os.system(cmd) 
    
    print "***********************************"
    print "PASS: Merge Files  ... DONE"
    print "***********************************"

