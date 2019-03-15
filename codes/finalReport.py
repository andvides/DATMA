#!/usr/bin/python
import os,argparse
import glob
import time   
        
def joinHtml():
    text='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"\n'
    text+='<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">\n'
    text+='<head>\n\t<meta charset="utf-8"/>\n</head>\n'
    text+='<body>\n'
    text+='<section>\n\t'
    text+='<h1>Bins report</h1>\n\t'
    text+='<iframe src="report.html" width="100%" frameborder="0" height="20%" ></iframe>\n'
    text+='</section>\n'
    text+='</section>\n'
    text+='<h1>Blastn report</h1>\n\t'
    text+='<iframe src="binsBlastn.html" width="100%" frameborder="0" height="600"></iframe>\n'
    text+='</section>\n'
    text+='</section>\n\t'
    text+='<h1>Kaiju report</h1>\n\t'
    text+='<iframe src="binsKaiju.html" width="100%" frameborder="0" height="600"></iframe>\n'
    text+="</section>\n</body>\n</html>\n"
    return text

def makeReport(typeReads,inputName):
    char='>'
    if (typeReads=='fastq'):
        char='@'
    count=0
    bases=0
    enable=True
    
    if (os.path.exists(inputName)):
        fsrc1 = open(inputName,'r') 
        for line in fsrc1:
            if enable:
                if(line.startswith(char)): #fasta or fastq
                    count+=1
                elif(line.startswith('+')): #fastq
                    enable=False
                else:
                    words=line.split("\n") 
                    words=words[0] 
                    bases+=len(words)
            else: #quality
                enable=True
        fsrc1.close()
        return str(count),str(bases)
    else:
        return str('Error'),str(inputName)
    
def tablaReport(section):
    if (section=='header'):
        text="<!DOCTYPE html>\n<html>\n<head>\n<style>\ntable, th, td {\n\t border: 1px solid black;\n\tborder-collapse: collapse;\n}\n</style>\n</head>\n<body>\n"
    elif (section=='title'):
        text="<table style='width:70%'>\n<tr>\n\t<th>Bin</th>\n\t<th>Size(reads)</th>\n\t<th>bp</th>\n\t<th>Contigs</th>\n\t<th>Genome</th>\n\t<th>ORFS</th>\n\t<th>bp</th>\n</tr>\n"
    elif(section=='end'):
        text="</table>\n</body>\n</html>\n"
    return text

# Read configuration file
def readConfigFile(fileName, param):
    error=True
    fsrc1 = open(fileName,'r') #abre el archivo 
    for line in fsrc1:
        words=line.split("\n") 
        words=words[0] 
        words=words.split()
        #basic parameters
        if line.startswith('-inputFile'): #mandatory!!!
            error=False
            inputFile=words[1]
            param.inputFile=inputFile
        elif line.startswith('-outputDir'):
            directory=words[1]
            param.outputFile=directory
        elif line.startswith('-cpus'):
            cpus=words[1]
            param.cpus=cpus
        elif line.startswith('-typeReads'): 
            typeReads=words[1]
            param.typeReads=typeReads 
            if not (param.typeReads=='sff' or param.typeReads=='illumina' or param.typeReads=='fastq' or param.typeReads=='fasta'):
                print "ERROR IN THE CONFIGURATION FILE"
                print "unsupported reads' format"
        elif line.startswith('-combine'): 
            param.combine='-c' 

def runCheckM (cpus,binDir,checkmDir):
    cmd='checkm lineage_wf -t '+str(cpus)+' '+binDir+' '+checkmDir
    print cmd
    os.system(cmd)
            
class parameters:
    manual = 'save the basic inputs'
    totalStages=8
    inputFile= ''
    outputFile='output'
    cpus=1
    typeReads='sff'      
    combine=' '


if __name__ == "__main__":
    #Final report
    description='DATMA: Distributed AuTomatic Metagenomc Assembly and Annotation framework report script \n'
    epilog='Authors: '
    epilog+='Benavides A, Alzate JF and Cabarcas F \n'
    
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-f', '--file', help='configuration file', required=True)

    args = parser.parse_args()
    param=parameters()
    
    if args.file :
        readConfigFile(args.file,param)    
        directory=param.outputFile
        
        #CheckM
        binDir=directory+'/bins/'
        checkmDir=directory+'/checkm_out/'
        runCheckM (param.cpus,binDir,checkmDir)
        
        #Run Krona
        combine=param.combine
        reportFile= directory+'/report.html'
        
        file = open(reportFile, 'w')
        file.write(tablaReport('header'))
        file.write(tablaReport('title'))
        file.close()
        
        pathDataset=directory+'/round_*/clameBin_*.f*'
        #print pathDataset
        finalBins=glob.glob(pathDataset)
        suffix=finalBins[0].split('.')
        typeReads=suffix[1]

        my_dict = {}
       
        pathDataset=directory+'/round_*/clameBin_*.'+typeReads
        #print pathDataset
        finalBins=glob.glob(pathDataset)
        for files in finalBins:
            suffix=files.split('.')
            suffix=suffix[0].split('/')
            suffix=suffix[-2]+'_'+suffix[-1]
            bins,bases=makeReport(typeReads,files)
            text=bins+","+bases+','
            my_dict[suffix]=text;
            #print suffix,':',text

        pathDataset=directory+'/bins/*contigs.fna'
        #print pathDataset
        finalBins=glob.glob(pathDataset)
        for files in finalBins:
            suffix=files.split('_contigs')
            suffix=suffix[0].split('/')
            suffix=suffix[-1]
            contigs,genome=makeReport('fasta',files)
            text=contigs+","+genome+','
            my_dict[suffix]+=text;
            #print suffix,'->',text
            
        pathDataset=directory+'/bins/*orfsFile.faa'
        print pathDataset
        finalBins=glob.glob(pathDataset)
        for files in finalBins:
            suffix=files.split('_orfsFile')
            suffix=suffix[0].split('/')
            suffix=suffix[-1]
            orfs,bases=makeReport('fasta',files)
            text=orfs+","+bases
            my_dict[suffix]+=text;
        
        #16S report
        suffix='all_16S'

        files=directory+'/16sSeq/all_16S.'+typeReads
        bins,bases=makeReport(typeReads,files)
        text=bins+","+bases
        my_dict[suffix]=text;
        
        text=",NA,NA,NA,NA"
        my_dict[suffix]+=text;
        
        
        #print '******************************************************'
        #print my_dict
        #print '******************************************************'
        sortednames=sorted(my_dict.keys(), key=lambda x:x.lower())
        print '******************************************************'
        print sortednames
        print '******************************************************'
        
        file = open(reportFile, 'a')
        for items in sortednames:
            #print '<tr>\n'
            file.write('<tr>\n')
            file.write('\t<th>'+items+'</th>\n')
            #print '\t<th>'+items+'</th>\n'
            for elements in my_dict[items].split(','):
                #print '\t<th>'+elements+'</th>\n'
                file.write('\t<th>'+elements+'</th>\n')
            file.write('<tr>\n')
            #print '<tr>\n'

        #file = open(reportFile, 'a')
        file.write(tablaReport('end'))
        file.close()
        
        pathDataset=directory+'/bins/*blastFile.tab'
        #print pathDataset
        finalBins=glob.glob(pathDataset)
        finalBins.sort()
        finalTaxo=' '#directory+'/16sSeq/contigs/16S_blastFile.tab '
        for taxos in finalBins:
            finalTaxo+=taxos+' '
        #print finalTaxo
        
        cmd= 'ktImportBLAST '+finalTaxo+' -o '+directory+'/binsBlastn.html '+combine
        print cmd
        os.system(cmd)

        pathDataset=directory+'/bins/*kaijuFile.txt'
        #print pathDataset
        finalBins=glob.glob(pathDataset)
        finalBins.sort()
        finalTaxo=' '#directory+'/16sSeq/contigs/16S_kaijuFile.krona '
        for taxos in finalBins:
            finalTaxo+=taxos+' '
        #print finalTaxo

        cmd= 'ktImportText '+finalTaxo+' -o '+directory+'/binsKaiju.html '+combine
        print cmd
        os.system(cmd)
        
        #16sReport
        cmd='ktImportRDP '+directory+'/16sSeq/all_16S.rdp'+' -o '+directory+'/16S_rdp.html '+combine
        print cmd
        os.system(cmd) 
    
        
        resumenFile= directory+'/resumenFile.html'
        file = open(resumenFile, 'w')
        file.write(joinHtml())
        file.close()
