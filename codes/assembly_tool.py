#!/usr/bin/python

from pycompss.api.task import task
from pycompss.api.parameter import *
"""
@task(binName=FILE, newblerOut=FILE_OUT, contigs=FILE_OUT, blastout=FILE_OUT, taxo=FILE_OUT)
def assembly(cpus,suffix,binName,newblerOut,contigs,blastout,taxo):
    import os
    import time
    import random
    a=random.randint(20, 100)
    """
    from shutil import copyfile
    copyfile(binName, newblerOut)
    copyfile(binName, contigs)
    copyfile(binName, blastout)
    copyfile(binName, taxo)
    """
    fos = open(newblerOut, 'w')
    fos.write('hi')
    fos.close()
    
    fos = open(contigs, 'w')
    fos.write('hi')
    fos.close()
    
    fos = open(blastout, 'w')
    fos.write('hi')
    fos.close()
    
    fos = open(taxo, 'w')
    fos.write('hi')
    fos.close()
    time.sleep( a )

"""
@task(binName=FILE, newblerOut=FILE_OUT, contigs=FILE_OUT, blastout=FILE_OUT, taxo=FILE_OUT)
def assembly(cpus,suffix,binName,newblerOut,contigs,blastout,taxo):
    import os
    from shutil import copyfile

    #5.3 using newbler to assemble the reads
    ensambleFolder='~/.Ensamblaje_'+suffix
    cmd="runAssembly -cpu "+cpus+" -o "+ensambleFolder+" "+binName+" > "+newblerOut
    if os.system(cmd):
        print  'program failed!->'+cmd
    else:
        #5.4 annotate
        contigfile=ensambleFolder+"/454LargeContigs.fna"
        xmlfile=ensambleFolder+"/454LargeContigs.xml"
        taxofile=ensambleFolder+"/454LargeContigs.taxo"
        annotate(cpus,contigfile,xmlfile,taxofile)
   
        #copy to the final files
        copyfile(contigfile, contigs)
        #cmd='cat '+contigfile+' > '+contigs
        #os.system(cmd)
        
        copyfile(xmlfile, blastout)
        #cmd='cat '+xmlfile+' > '+blastout
        #os.system(cmd)
        
        copyfile(taxofile, taxo)
        #cmd='cat '+taxofile+' > '+taxo
        #os.system(cmd)
        
        #remove temp files
        cmd='rm -r '+ensambleFolder+'*'
        os.system(cmd)

#6.Annotation for each bin
def annotate(cpus,contigfile,xmlfile,taxofile):
    import os
    
    # Read database
    #database='/home/users/andresb/database.txt'
    #db = open(database,'r').read().split('\n')
    #cmd='blastn -query '+contigfile+' -out '+xmlfile+' -db  '+str(db[0])+' -outfmt 5 -num_alignments 5 -num_threads '+cpus
    
    #blast
    cmd='blastn -query '+contigfile+' -out '+xmlfile+' -db /home/db/NT/NT_May_2017/nt -outfmt 5 -num_alignments 5 -num_threads '+cpus
    print cmd
    os.system(cmd)

    #xml2tab
    cmd='xml2tabV2.py'+" -xml "+xmlfile+" > "+taxofile
    print cmd
    os.system(cmd)
