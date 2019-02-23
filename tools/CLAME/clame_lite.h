/*---------------------------------------------------------------
 *
 *  CLAME:"Clasificador Metagenomico"
 *	
 *  CLAME is a binning software for metagenomic reads.
 *  It immplements a fm-index search algorithm for nucleotide 
 *  sequence alignment. Then it uses strongly connected component strategy
 *  to bin the similar sequences.
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
 *  Version 2.2
 *  Binning method was modified to include edges parameter with non inclusion strategy 
 *  
 ---------------------------------------------------------------*/

#include <sdsl/suffix_arrays.hpp>
#include "iostream"
#include<fstream>
#include <ctime>

#include <fcntl.h>
#include <sstream> 
#include <vector>
#include <omp.h>
#include <string>
#include <algorithm>
#include <iomanip>

#include<cstdint>
#include <cstdlib>
#include <tr1/unordered_map>
#include <vector>


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
//using std::unordered_map;
using std::tr1::unordered_map;

using namespace sdsl;
using namespace std;

struct Args {bool multiFasta; bool fastq; bool outputFile; bool numT; bool bases_Threshold; bool print; bool fm9; bool edges; bool tol; bool sizeBin; bool ld; bool w; bool forward_reverse;};
struct Names {string multiFasta; string outputFile; string fm9;};
struct Parameters {float edges; float tol; int numThreads; int query_size; bool enablePrint; bool loadFM9;int sizeBin;bool fastq; int  ld; int w; string forward_reverse;};

string reverse(string str);
int mstrlen(const char arg[]);
void printerror(const char arg[]);
bool IsParam(char arg[],const char comp[]);
bool readArguments(int argc,char *argv[], Args *args,Names *names, Parameters *parameters);
bool readFasta(Names *names,vector<string> *bases,string *fasta,vector<uint32_t> *index, vector<string> *title);
bool readFastQFile(Names *names,vector<string>* bases, string* fasta, vector<uint32_t> *index, vector<string> *title, vector<string> *qual);
bool alignemnt(Names *names,Parameters *parameters,vector<string> *bases,string *fasta, vector<uint32_t> *index, int* queryList, vector<int>* MatrixList);
void printResult(Names *names, int numberOFreads, vector<int> *MatrixList);

uint32_t indx2Loc(uint32_t locs, vector<uint32_t> *index);
void binning(Names *names, Parameters *parameters, vector<string> *title, int *queryList, vector<int> *MatrixList, int numberOFreads);
void binningPrint(Names *names, Parameters *parameters, vector<string> *title, int *queryList, vector<int> *MatrixList, int numberOFreads, vector<string> *bases,vector<string> *qual);
void static_parameters(int *stack, int size, vector<int> *MatrixList, int *queryList, int numBin, int *maximun, int *minimum, float cut, float *stadistic) ;
