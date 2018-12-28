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

// I N C L U D E S ***********************************************************/
#include <zlib.h>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include<fstream>
#include <ctime>
#include <fcntl.h>
#include <sstream> 
#include <vector>
#include <string.h>
#include <unistd.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/* D E F I N E S *************************************************************/
#define VERSION "V2"
#define PRG_NAME "rapifilt"
#define FASTQ_FILENAME_MAX_LENGTH 1024
#define input_FILENAME_MAX_LENGTH 1024
#define maxReadLen 2000
using namespace std;

// G L O B A L S *************************************************************/
struct Args 
{
    bool outputFile;    // Output file name 
    bool illuminaFile;  // Input illumina file name 
    bool sffFile;       // Input illumina file name  
    bool gz;            // Input illumina gz file name  
    bool fastq;         // Input simple fastq file name  
    bool enableFasta;   // selec between fasta(1) and fastq output(0)
    int left_cut;       // left cut value 
    int right_cut;      // right cut value 
    int windows_size;   // windows size to check 
    int min_len;        // Filter sequence shorter than min_len
    int max_len;        // Filter sequence larger than max_len
    int tb;             //removes n bases from the begins of sequencing fragments (int default 0)
    int te;             //removes n bases from the ends of sequencing fragments (int default 0)
    int binSize;        //to compute statistic per base (int default 1)
};

struct Names 
{
    string outputFile; 
    string gz_illuminaFile1; 
    string gz_illuminaFile2; 
    string illuminaFile1; 
    string illuminaFile2; 
    string sffFile; 
    string fastq;
    
};

/* P R O T O T Y P E S *******************************************************/
void process_options(int argc, char *argv[], struct Args * args, struct Names * names);
void help_message(void);
int mstrlen(const char arg[]);
bool IsParam(char arg[],const char comp[]);
void GC_content(char *bases1, int size, int * Gbase, int  * Cbase, int bouderLeft1, int bouderRight1);
void statistic_mean (unsigned int *quality1_int, int binSize,int size, int *max_read_len, int *maxQ, int *minQ, int *total_reads, int *qual_position, int bouderLeft1, int bouderRight1);
//void statistic_mean (vector<unsigned int> *quality1_int, int binSize,int size, int *max_read_len, int *maxQ, int *minQ, int *total_reads, int *qual_position, int bouderLeft1, int bouderRight1);
