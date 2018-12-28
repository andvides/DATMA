/*---------------------------------------------------------------
 *
 *  RAPIFILT:"RAPId FILTer"
 *	
 *  RAPIFILT is a fast computational tool designed to trim sequences 
 *  using the quality scores of bases within individual read.
 *  
 *  Copyright (C) 2017 Benavides A.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/ .
 *   
 *  Version 2.0
 ---------------------------------------------------------------*/

#include "generic.h"

void process_options(int argc, char *argv[], struct Args * args, struct Names * names)
{
    bool   error = false;
    int    argum = 1;
    if(argc==1)
        help_message(), exit(0);
    else
        while (argum < argc && not error) 
        {
            error = true;
            
            if ( IsParam(argv[argum],"-h") ) 
            {
                help_message();
                exit(0);
            }
            
            if ( IsParam(argv[argum],"-v") ) 
            {
                cout <<PRG_NAME<<"_"<<VERSION<<endl;
                exit(0);
            }
            
            if ( IsParam(argv[argum],"-o") ) //output file
            {
                error = false;
                argum ++;
                if (argum < argc)  
                {
                    args->outputFile = true;
                    names->outputFile = argv[argum];
                } 
                else 
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-i") ) //input Illumina file
            {
                  
                argum ++;
                if (argum < argc)  
                {
                    names->illuminaFile1 = argv[argum];
                    argum ++;
                    if (argum < argc)
                    {
                        error = false;
                        args->illuminaFile = true;
                        names->illuminaFile2 = argv[argum];
                    }
                } 
                else 
                    error = true;
                    argum ++;
                continue;
            }
            
            
            if ( IsParam(argv[argum],"-gz") ) //input gz Illumina file
            {
                argum ++;
                if (argum < argc)  
                {
                    names->gz_illuminaFile1 = argv[argum];
                    argum ++;
                    if (argum < argc)
                    {
                        error = false;
                        args->gz= true;
                        names->gz_illuminaFile2 = argv[argum];
                    }
                } 
                else 
                    error = true;
                argum ++;
                continue;
            }   
            
            if ( IsParam(argv[argum],"-fastq") ) //input 454 file
            {
                argum ++;
                if (argum < argc)  
                {
                    error = false;
                    args->fastq = true;
                    names->fastq = argv[argum];
                } 
                else 
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-sff") ) //input 454 file
            {
                argum ++;
                if (argum < argc)  
                {
                    error = false;
                    args->sffFile = true;
                    names->sffFile = argv[argum];
                } 
                else 
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-f") ) 
            {
                error = false;
                argum ++;
                args->enableFasta = true;
                continue;
            }
            
            if ( IsParam(argv[argum],"-l") ) 
            {
                error = false;
                argum ++;
                if (argum < argc)  
                    args->left_cut = atoi(argv[argum]);
                else 
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-r") ) 
            {
                error = false;
                argum ++;
                if (argum < argc)  
                    args->right_cut = atoi(argv[argum]);
                else 
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-w") ) 
            {
                error = false;
                argum ++;
                if (argum < argc)  
                    args->windows_size = atoi(argv[argum]);
                else 
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-m") ) 
            {
                error = false;
                argum ++;
                if (argum < argc)  
                    args->min_len = atoi(argv[argum]);
                else 
                    error = true;
                argum ++;
                continue;
            }

            if ( IsParam(argv[argum],"-mx") )
            {
                error = false;
                argum ++;
                if (argum < argc)
                    args->max_len = atoi(argv[argum]);
                else
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-tb") ) 
            {
                error = false;
                argum ++;
                if (argum < argc)  
                    args->tb = atoi(argv[argum]);
                else 
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-te") ) 
            {
                error = false;
                argum ++;
                if (argum < argc)  
                    args->te = atoi(argv[argum]);
                else 
                    error = true;
                argum ++;
                continue;
            }
            
            if ( IsParam(argv[argum],"-bin") ) 
            {
                error = false;
                argum ++;
                if (argum < argc)  
                    args->binSize = atoi(argv[argum]);
                else 
                    error = true;
                argum ++;
                continue;
            }
    }

    if (error)
        help_message(), exit(0);
}

void help_message() 
{
    cout << "RAPIFILT:RAPId FILTer" << endl;
    cout << "version 2.2 Decemebr 2017"<<endl;
    cout << "Authors"<<endl;
    cout << "Benavides A, Alzate JF and Cabarcas F"<<endl;
    cout<< "Usage: "<<PRG_NAME<<" [options]"<<endl;
    
    cout<< "-h"<<"\t\t\tThis help message\n";
    cout<< "-v"<<"\t\t\tProgram and version information\n";
    cout<< "-f"<<"\t\t\tEnable fasta output (default fastq)\n";
    cout<< "-l<int>"<<"\t\t\tSet lef-cut value for quality scores (int default 0)\n";
    cout<< "-r<int>"<<"\t\t\tSet right-cut value for quality scores (int default 0)\n";
    cout<< "-w<int>"<<"\t\t\tSet windows size to check on the quality scores (int default 1)\n";
    cout<< "-m<int>"<<"\t\t\tFilter sequence shorter than min_len (int default 1)\n";
    cout<< "-mx<int>"<<"\t\tFilter sequence larger than max_len (int default 5000)\n";
    cout<< "-fastq<fastq file>"<<"\tsingle fastq input (file.fastq)\n";
    cout<< "-sff<454 files>"<<"\t\t454 input (file.sff)\n";
    cout<< "-i<illumina files>"<<"\tIllumina inputs(file1.fastq file2.fastq)\n";
    cout<< "-gz<gzipped file>"<<"\tgzipped file(file1.fastq.gz file2.fastq.gz)\n";
    cout<< "-o<fastq_file>"<<"\t\tDesired fastq output file."<<" If not specified to stdout\n";
    cout<< "-tb<int>"<<"\t\tRemove n bases from the begins of sequencing fragments (int default 0)\n";
    cout<< "-te<int>"<<"\t\tRemove n bases from the ends of sequencing fragments (int default 0)\n";
    cout<< "-bin<int>"<<"\t\tBin size used to compute statistic per base (int default 1)\n";
    cout<< "Quality input and output is calculated according the formula taken from http://maq.sourceforge.net/fastq.shtml"<<endl;	
}

void GC_content(char *bases1, int size, int * Gbase, int  * Cbase, int bouderLeft1, int bouderRight1)
{
    //GC-content
    for (int i=bouderLeft1; i<size-bouderRight1; i++)
        if(bases1[i]=='c' || bases1[i]=='C')
            *Cbase++;
        else if(bases1[i]=='g' || bases1[i]=='G')
            *Gbase++;                       
}

void statistic_mean (unsigned int *quality1_int, int binSize,int size, int *max_read_len, int *maxQ, int *minQ, int *total_reads, int *qual_position, int bouderLeft1, int bouderRight1)
{
    if(size>*max_read_len)
        *max_read_len=size;
    
    //Qualiti score across all bases 
    int i=0, index=0;
    for (i=bouderLeft1; i<size; i++)
    {
        index=i/binSize;
        maxQ[index]=(maxQ[index]>quality1_int[i])?maxQ[index]:quality1_int[i];
        minQ[index]=(minQ[index]<quality1_int[i])?minQ[index]:quality1_int[i];
        total_reads[index]+=1;
        qual_position[index]+=quality1_int[i];
    }
}
/*void statistic_mean (vector<unsigned int> *quality1_int, int binSize,int size, int *max_read_len, int *maxQ, int *minQ, int *total_reads, int *qual_position, int bouderLeft1, int bouderRight1)
{
    if(size>*max_read_len)
        *max_read_len=size;
    
    //Qualiti score across all bases 
    int i=0, index=0;
    for (std::vector<unsigned int>::iterator it = quality1_int->begin()+bouderLeft1; it != quality1_int->end()-bouderRight1; ++it)
        index=i/binSize,maxQ[index]=(maxQ[index]>*it)?maxQ[index]:*it,minQ[index]=(minQ[index]<*it)?minQ[index]:*it,total_reads[index]+=1, qual_position[index]+=*it, i++;
                    
                    
}*/

int mstrlen(const char arg[]) 
{
    int out = 0;
    while(arg[out] != 0)
        out++;
    return out;
}

bool IsParam(char arg[],const char comp[])
{
    int len_comp = mstrlen(comp);
    bool equal = true;
    if (mstrlen(arg) == len_comp) 
    {
        for (int i = 0; i < len_comp; i++) 
            if (arg[i] != comp[i]) 
                equal = false;
    } 
    else 
        equal = false;
  
    return equal;
}
