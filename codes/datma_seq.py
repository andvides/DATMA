#!/usr/bin/python
#
#runcompss --summary --master_name=aletoso -d --project=./project.xml --resources=./resources.xml --lang=python --master_port=43100 ./clameWorkFlow.py -f simple.txt
#
import os,argparse
import glob
import time

#CLAMEWORFLOW MODULES
from general import *
from removal_tool import seq_16Sremoval
from filterQuality import cleanReads
from binning_tool import *
from assembly_annotation_tool_seq import *
#from pycompss.api.api import waitForAllTasks

  
if __name__ == "__main__":

    description='DATMA: Distributed AuTomatic Metagenomc Assembly and Annotation framework \n'
    description+='Version 3.0 December 2018 \n'
    epilog='Authors: '
    epilog+='Benavides A, Sanchez F, Alzate JF and Cabarcas F \n'
    
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-f', '--file', help='configuration file', required=True)
    args = parser.parse_args()
        
    if args.file:
        start_time = time.time()

        fileName=args.file

        param=parameters()
        tool=tools()
        outName=outputNames()
        direct=directories()
        flags=stageFlags()
        
        #Check for all the tools
        is_tool(tool)
    
        #Read configuration file
        readConfigFile(fileName, param)    
        
        #Update pipeline stages execution
        pipelineFlow(param,flags)
        
        #build output directory
        builDirectory(direct,param,flags)
        
        #1.Clean de reads 
        cleanReads(direct,param,outName,flags)
        
        #2.Remove 16S sequences
        seq_16Sremoval(direct,param,outName,flags)
        
        #CLAME rounds
        i=0
        dirBase=direct.rounds
        bases=param.bases.split(',')
        
        

        for items in bases:
            param.base=bases[i]
            print 'Current round ',i,'->',param.base
            
            #make direct for the current round
            direct.rounds=dirBase+str(i)+"_b"+str(param.base)
            directory=direct.globalOutput+'/'+direct.rounds
            
            cmd="mkdir "+directory
            print cmd
            os.system(cmd)
            
            #update the file to round
            file2bin=direct.globalOutput+'/'+outName.balance+'.'+param.typeReads
            cmd="mv "+file2bin+" "+directory
            print cmd
            os.system(cmd)
            
            
            #3.Binning the clean reads <CLAME>
            BinReads(direct,param,outName,flags,i)
            i+=1
           
            
        #4. Assemble and Annotate each bin by separated
        direct.rounds=dirBase
        directory=direct.globalOutput+'/'+direct.rounds
        pathDataset = directory+'*/'+outputNames.clame+'_*.fast*'  
        bins=glob.glob(pathDataset)
        bins.sort(key=lambda x: int(x.split(outputNames.clame)[1].split(".")[0][1:]))
        print pathDataset
        print bins
        binNum=0
        
        binDir=direct.globalOutput+'/bins/'
        cmd="mkdir "+binDir
        print cmd
        os.system(cmd)
            
        for binName in bins: #for each bin generated
            
            suffix=binName.split('.')
            prefix=suffix[1]
            suffix=suffix[0].split('/')
            suffix=suffix[-2]+'_'+suffix[-1]
            
            contigsFile=binDir+suffix+'_contigs.fna'
            orfsFile=binDir+suffix+'_orfsFile.faa'
            blastFile=binDir+suffix+'_blastFile.tab'
            database_nt=param.database_nt
            kaijuFile=binDir+suffix+'_kaijuFile.txt'
            database_kaiju=param.database_kaiju
            reportFileTemp=binDir+suffix+'_report.txt'
            
            assembly_annotation(str(suffix),str(prefix),fileName,binName,contigsFile,orfsFile,blastFile,str(database_nt),kaijuFile,str(database_kaiju),str(param.cpus),param.assembly)
    
        print 'END DATMA WORK FLOW'
        print time.time() - start_time, "seconds"
        
    else:
        print description
        print epilog
        
        

