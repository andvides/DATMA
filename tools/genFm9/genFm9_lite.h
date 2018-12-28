/*---------------------------------------------------------------
 *
 *   CLAME:"Clasificador Metagenómico"
 *
 *   CLAME is a binning software for metagenomic reads.
 *   It immplements a fm-index search algorithm for nucleotide 
 *   sequence alignment. Then it uses strongly connected component strategy
 *   to bin the similar sequences.
 *   
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

#include <sys/resource.h>


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
//using std::unordered_map;
using std::tr1::unordered_map;

using namespace sdsl;
using namespace std;

struct Args {bool multiFasta; bool fastq; bool outputFile; bool numT; bool bases_Threshold; bool print; bool fm9; bool lu; bool ld; bool sizeBin;};
struct Names {string multiFasta; string outputFile; string fm9;};
struct Parameters {int ld; int lu; int numThreads; int query_size; bool enablePrint; bool loadFM9;int sizeBin;};

string reverse(string str);
int mstrlen(const char arg[]);
void printerror(const char arg[]);
bool IsParam(char arg[],const char comp[]);
bool readArguments(int argc,char *argv[], Args *args,Names *names, Parameters *parameters);
bool readFasta(Names *names,vector<string> *bases,string *fasta,vector<uint32_t> *index,vector<string> *title);
bool readFastQFile(Names *names,vector<string>* bases, string* fasta, vector<uint32_t> *index,vector<string> *title);
bool alignemnt(Names *names,Parameters *parameters,vector<string> *bases,string *fasta, vector<uint32_t> *index, int* queryList, vector<int>* MatrixList);
void printResult(Names *names, int numberOFreads, vector<int> *MatrixList);

uint32_t indx2Loc(uint32_t locs, vector<uint32_t> *index);
void binning(Names *names, Parameters *parameters, vector<string> *title, int *queryList, vector<int> *MatrixList, int numberOFreads);
void binningPrint(Names *names, Parameters *parameters, vector<string> *title, int *queryList, vector<int> *MatrixList, int numberOFreads, vector<string> *bases);

