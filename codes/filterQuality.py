#!/usr/bin/python
import os,sys
from merge_tool import *

def cleanReads(direct,param,outName,flags):
    #Input-ouput names
    directory= direct.globalOutput+'/'+direct.filterQuality
    inputFile=param.inputFile
    outputFile=directory+'/'+outName.clean

    cpus=param.cpus
    typeReads=param.typeReads
    cleanTool=str(param.cleanTool)
    trimmomatic_path=param.trimmomatic_path

    te=str(param.te)
    tb=str(param.tb)
    lq=str(param.lq)
    rq=str(param.rq)
    m=str(param.m)
    w=str(param.wq)
    aux_param=str(param.quality_aux)

    if flags.quality:
        #make directory to save the results
        cmd="mkdir "+directory
        print(cmd)
        os.system(cmd)
        
        #Quality filter
        if typeReads=='sff':
            if(cleanTool=='rapifilt'):
                cmd='rapifilt -sff '+inputFile+' -o '+outputFile+' -tb '+tb+' -te '+te+' -l '+lq+' -r '+rq+' -m '+m+' -w '+w+' '+aux_param
                print(cmd)
                os.system(cmd) 
                #update the format for the next stages
                param.typeReads='fastq'
            else:
                print("ERROR: Format not supported by the quality control tool")
                sys.exit(10)
            cmd='fastqc '+directory+'/'+outName.clean+'.fastq'+' -o '+directory+'/'
            print(cmd)
            os.system(cmd) 
        elif typeReads=='illumina':
            files=inputFile.split(',')
            #update the format for the next stages
            param.typeReads='fastq'
            if(cleanTool=='prinseq'):
                cmd='prinseq-lite.pl -fastq  '+files[0]+' -fastq2 '+files[1]+' -out_format 3 -no_qual_header -out_good '+outputFile+' -out_bad null -min_len '+m+' -rm_header -trim_qual_right '+rq+' -trim_qual_left '+lq+' -trim_qual_type min  -trim_qual_window '+w+' '+aux_param
                print(cmd)
                os.system(cmd) 
            elif(cleanTool=='trimmomatic'):
                if(trimmomatic_path=='none'):
                    cmd='trimmomatic'+' PE -phred33 -threads '+cpus+' '+files[0]+' '+files[1]+' '+outputFile+'.'+param.typeReads+' LEADING:'+lq+' TRAILING:'+rq+'  MINLEN:'+m+' '+aux_param
                else:
                    cmd='java -jar '+trimmomatic_path+'trimmomatic'+' PE -phred33 -threads '+cpus+' '+files[0]+' '+files[1]+' '+outputFile+'.'+param.typeReads+' LEADING:'+lq+' TRAILING:'+rq+'  MINLEN:'+m+' '+aux_param
                print(cmd)
                os.system(cmd) 
            elif(cleanTool=='fastx'): #Fastx: NO pair end ouput, no left cut, no windows otion Q33 format
                print("ERROR: Format not supported by the quality control tool")
                sys.exit(10)
            else : #(cleanTool=='rapifilt') default tool    
                cmd='rapifilt -i '+files[0]+' '+files[1]+' -o '+outputFile+' -tb '+tb+' -te '+te+' -l '+lq+' -r '+rq+' -m '+m+' -w '+w+' '+aux_param
                print(cmd)
                os.system(cmd) 
                
            #merge files
            temName='mergeReads'
            merge(param,outputFile,temName,directory)
            #rename the temp files to the regular names
            cmd='mv '+directory+'/'+temName+'.extendedFrags.fastq '+directory+'/'+outName.clean+'.fastq'
            print(cmd)
            os.system(cmd)
            
            cmd='fastqc '+directory+'/'+outName.clean+'.fastq'+' -o '+directory+'/'
            print(cmd)
            os.system(cmd) 
        elif typeReads=='fastq':
            if(cleanTool=='prinseq'):
                cmd='prinseq-lite.pl -fastq  '+inputFile+' -out_format 3 -no_qual_header -out_good '+outputFile+' -out_bad null -min_len '+m+' -rm_header -trim_qual_right '+rq+' -trim_qual_left '+lq+' -trim_qual_type min  -trim_qual_window '+w+' '+aux_param
                print(cmd)
                os.system(cmd) 
            elif(cleanTool=='fastx'):#Fastx: NO pair end ouput, no left cut, no windows otion Q33 format
                cmd='fastq_quality_trimmer -t '+rq+' -i '+inputFile+' -o '+outputFile+'.fastq -Q 33 -l '+m+' '+aux_param
                print(cmd)
                os.system(cmd)  
            elif(cleanTool=='trimmomatic'):
                if(trimmomatic_path=='none'):
                    cmd='trimmomatic'+' SE -phred33 -threads '+cpus+' '+inputFile+' '+outputFile+'.'+param.typeReads+' LEADING:'+lq+' TRAILING:'+rq+'  MINLEN:'+m+' '+aux_param
                else:
                    cmd='java -jar '+trimmomatic_path+'trimmomatic'+' SE -phred33 -threads '+cpus+' '+inputFile+' '+outputFile+'.'+param.typeReads+' LEADING:'+lq+' TRAILING:'+rq+'  MINLEN:'+m+' '+aux_param
                print(cmd)
                os.system(cmd) 
            else:#(cleanTool=='rapifilt') default tool    
                cmd='rapifilt -fastq '+inputFile+' -o '+outputFile+' -tb '+tb+' -te '+te+' -l '+lq+' -r '+rq+' -m '+m+' -w '+w+' '+aux_param
                print(cmd)
                os.system(cmd)  
        
            cmd='fastqc '+directory+'/'+outName.clean+'.fastq'+' -o '+directory+'/'
            print(cmd)
            os.system(cmd) 
        elif typeReads=='fasta':
            #symbolic link to the original file
            os.symlink(inputFile,outputFile+'.fasta')
            print('NO quality file type')
        else:
            print('Quality Filter stage:ERROR IN THE CONFIGURATION FILE')
            print("Unknown reads' format")
            sys.exit(0)
        
        print("***********************************")
        print("PASS: Quality filter stage ... DONE")
        print("***********************************")
    else:
        outputFile+='.'+typeReads
        if not (os.path.exists(outputFile)):
            if typeReads=='sff':
                print('Quality Filter stage:ERROR IN THE CONFIGURATION FILE')
                print("Filter stage madatory for this format") 
                sys.exit(10)
            elif typeReads=='illumina':
                files=inputFile.split(',')
                param.typeReads='fastq'    #update the format for the next stages
                outputFile=directory+'/'+outName.clean+'.fastq'
                if not (os.path.exists(outputFile)):
                    #merge files
                    temName='mergeReads'
                    mergeGeneric(param,files[0],files[1],temName,directory)
                    #rename the temp files to the regular names
                    cmd='mv '+directory+'/'+temName+'.extendedFrags.fastq '+outputFile
                    print(cmd)
                    os.system(cmd)
            else:
                #make directory to save the results
                cmd="mkdir "+directory
                print(cmd)
                os.system(cmd)
                #symbolic link to the original file
                os.symlink(inputFile,outputFile)
        
        print("***********************************")
        print("PASS: NO Quality filter stage      ")
        print("***********************************")
