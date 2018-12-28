#!/usr/bin/python

from pycompss.api.task import task
from pycompss.api.parameter import *

@task(binName=FILE, newblerOut=FILE_OUT, contigs=FILE_OUT, blastout=FILE_OUT, taxo=FILE_OUT)
def annotate(cpus,suffix,binName,newblerOut,contigs,blastout,taxo):
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



