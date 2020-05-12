#!/usr/bin/python
#
#DATMA: Distributed AuTomatic Metagenomc Assembly and Annotation framework          #
#Version 4.0 November 2019                                                          #
#Authors:                                                                           #
#Benavides A, Sanchez F, Alzate JF, and Cabarcas F                                  #
#
import os,argparse, errno
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
    description+='Version 4.0 November 2019 \n'
    epilog='Authors: '
    epilog+='Benavides A, Sanchez F, Alzate JF, and Cabarcas F \n'
    
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
        
        startStage=param.startStage
        print("startStage",startStage)
        if(startStage>4):
            print("END first part")
            exit(10)

        #Update pipeline stages execution
        pipelineFlow(param,flags)
        
        #build output directory
        builDirectory(direct,param,flags)
        
        #1.Clean de reads 
        illuminaReads=param.typeReads;
        if(param.endWorkflow==1):
            print("End in Quality control")
            exit(10)
        cleanReads(direct,param,outName,flags)
        print("Quality control")
        print(time.time() - start_time, "seconds")
    
        #2.Remove 16S sequences
        if(param.endWorkflow==2):
            print("End in 16S ribosomal")
            exit(10)
        seq_16Sremoval(direct,param,outName,flags)
        print("16S  ribosomal")
        print(time.time() - start_time, "seconds")
        
        #3. CLAME rounds
        if(param.endWorkflow==3):
            print("End in Binning stage")
            exit(10)
        
        dirBase=direct.rounds
        if flags.clame: #start in binning stage
            i=0
            bases=param.bases.split(',')
            
            #Removing a possible exclude list
            list2Exclude=direct.globalOutput+'/'+outName.balance+'.exclude'
            cmd='rm -f '+list2Exclude
            print(cmd)
            os.system(cmd) 
            
            #making or using global FM9
            if(param.fm9=="unassigned"):
                file2bin=direct.globalOutput+'/'+outName.balance+'.'+param.typeReads
                outputFile=direct.globalOutput+'/'+outName.balance
                genFM9(param.cpus,file2bin,outputFile,param.typeReads,param.block)
            else:
                outputFile=direct.globalOutput+'/'+outName.balance
                j=0
                for file in glob.glob(param.fm9+"*.fm9"):
                    target=param.fm9+str(j)+'.fm9'
                    link_name=outputFile+str(j)+'.fm9'
                    if not (target==outputFile+str(j)+'.fm9'):
                        print(target,"----->",link_name)
                        try:
                            os.symlink(target, link_name)
                        except OSError as e:#except OSError, e:
                            if e.errno == errno.EEXIST:
                                os.remove(link_name)
                                os.symlink(target, link_name)
                            else:
                                raise e
                    j+=1
                
            for items in bases:
                param.base=bases[i]
                print('Current round ',i,'->',param.base)
                
                #make direct for the current round
                direct.rounds=dirBase+str(i)+"_b"+str(param.base)
                directory=direct.globalOutput+'/'+direct.rounds
                
                cmd="mkdir -p "+directory
                print(cmd)
                os.system(cmd)
                        
                #3.Binning the clean reads <CLAME>
                BinReads(direct,param,outName,flags,i)
                i+=1
            print("Binning")
            print(time.time() - start_time, "seconds")
       
        #4. Assemble and Annotate each bin by separated
        if(param.endWorkflow==4 or startStage >4):
            print("End in assembly stage")
            exit(10)
            
        direct.rounds=dirBase
        directory=direct.globalOutput+'/'+direct.rounds
        #pathDataset = directory+'*/'+outputNames.clame+'_*.fast*'  
        pathDataset = directory+'*/'+outputNames.clame+'_*.list'  
        bins=glob.glob(pathDataset)
        bins.sort(key=lambda x: int(x.split(outputNames.clame)[1].split(".")[0][1:]))
        print(pathDataset)
        print(bins)
        binNum=0
        
        #if(param.endWorkflow==5):
            #print "End in annotation stage"
            #exit(10)
        binDir=direct.globalOutput+'/bins/'
        cmd="mkdir -p "+binDir
        print(cmd)
        os.system(cmd)
        
        #pass the list 2 rawReads
        param.typeReads=illuminaReads
        for binName in bins: #for each bin generated
            list2rawReads(binName,param,direct.globalOutput)
        print("End pass the list 2 rawReads")
        
        
        """pathDataset = directory+'*/'+outputNames.clame+'_*.fast*'  
        bins=glob.glob(pathDataset)
        bins.sort(key=lambda x: int(x.split(outputNames.clame)[1].split(".")[0][1:]))
        print(pathDataset)
        print(bins)
        binNum=0
        """
        
        for binName in bins: #for each bin generated

            fileName=binName.split('.')
            suffix=binName.split('.')
            prefix=param.typeReads#suffix[1]
            suffix=suffix[0].split('/')
            suffix=suffix[-2]+'_'+suffix[-1]
            
            assembler=param.assembly
            use_ref=param.use_ref;
            ref_name=param.ref_name;
            asm_aux=param.asm_aux
            ORF_aux=param.ORF_aux
            blast_aux=param.blast_aux
            kaiju_aux=param.kaiju_aux
            quast_aux=param.quast_aux
            
            contigsFile=binDir+suffix+'_contigs.fna'
            contigsFile_stat=binDir+suffix+'_contigs_qual.html'
            orfsFile=binDir+suffix+'_orfsFile.faa'
            blastFile=binDir+suffix+'_blastFile.tab'
            database_nt=param.database_nt
            kaijuFile=binDir+suffix+'_kaijuFile.txt'
            database_kaiju=param.database_kaiju
            reportFileTemp=binDir+suffix+'_report.txt'
            
            if prefix=='illumina':
                forwardFile=fileName[0]+'_1.fastq'
                reverseFile=fileName[0]+'_2.fastq'
            else:
                forwardFile=fileName[0]+'.'+prefix
                reverseFile=fileName[0]+'.'+prefix
                
            print("Files to tranfer:")
            print(str(suffix),str(prefix),forwardFile,reverseFile,contigsFile,contigsFile_stat,orfsFile,blastFile,str(database_nt),kaijuFile,str(database_kaiju),str(param.cpus),str(assembler),bool(use_ref),str(ref_name),str(asm_aux),str(ORF_aux),str(blast_aux),str(kaiju_aux),str(quast_aux))
            print("---------------------------------------")
            
            assembly_annotation(str(suffix),str(prefix),forwardFile,reverseFile,contigsFile,contigsFile_stat,orfsFile,blastFile,str(database_nt),kaijuFile,str(database_kaiju),str(param.cpus),str(assembler),bool(use_ref),str(ref_name),str(asm_aux),str(ORF_aux),str(blast_aux),str(kaiju_aux),str(quast_aux))
            
            
    
        print('END DATMA WORK FLOW')
        print(time.time() - start_time, "seconds")
        
    else:
        print(description)
        print(epilog)
