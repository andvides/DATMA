
/*  Select fastQ reads from list"
 *
 *        -h --help
 *        -fastq  FILE
 *        -lista  FILE
 *        -fasta  FILE
 *        -fastq2fasta
 *                    */

#include<iostream>
#include<cstdint>
#include<fstream>
#include<string>
#include <cstdlib>
#include <tr1/unordered_map>
#include <ctime>
#include <fcntl.h>
#include <string.h>
#include <sstream> 
#include <unistd.h>


using std::tr1::unordered_map;
using namespace std;

struct Args { bool fastq; bool lista; bool randn; bool fastq2fasta; bool fasta; bool fasta_sel;};
struct Names { string fastq; string lista; int randn; bool fastq2fasta; string fasta; bool fasta_sel;};

int mstrlen(const char arg[]);
void printerror(const char arg[]);
bool IsParam(char arg[],const char comp[]);
void selectfastq(Names& names);
void selectfasta(Names& names);
void fastq2fasta(Names& names);
void getrandom(Names& names);
static uintmax_t wc(char const *fname);

int main(int argc,char *argv[])
{
  bool   error = false, found1;
  Args   args = {false,false,false,false,false,false};
  Names  names;
  int    argum = 1;

  while (argum < argc && not error) {
    found1 = false;
    if ( IsParam(argv[argum],"-fastq") ) {
      argum ++;
      if (argum < argc) { //verifica si el parametro existe
        args.fastq = true;
        names.fastq = argv[argum];
      } else {
        error = true;
      }
      argum ++;
      found1 = true;
      continue;
    }
    if ( IsParam(argv[argum],"-list") ) {
      argum ++;
      if (argum < argc) { //verifica si el parametro existe
        args.lista = true;
        names.lista = argv[argum];
      } else {
        error = true;
      }
      argum ++;
      found1 = true;
      continue;
    }
    if ( IsParam(argv[argum],"-fasta") ) {
      argum ++;
      if (argum < argc) { //verifica si el parametro existe
        args.fasta = true;
        names.fasta = argv[argum];
      } else {
        error = true;
      }
      argum ++;
      found1 = true;
      continue;
    }
    if ( IsParam(argv[argum],"-fasta_sel") ) {
      args.fasta_sel = true;
      found1 = true;
      argum ++;
      continue;
    }
    if ( IsParam(argv[argum],"-fastq2fasta") ) {
      args.fastq2fasta = true;
      found1 = true;
      argum ++;
      continue;
    }
    if ( IsParam(argv[argum],"-h") ) {
      error = true;
      found1 = true;
      continue;
    }
    if ( IsParam(argv[argum],"-random") ) {
      argum ++;
      if (argum < argc) { //verifica si el parametro existe
        args.randn = true;
        names.randn = atoi(argv[argum]);
      } else {
        error = true;
      }
      argum ++;
      found1 = true;
      continue;
    }
    if (not found1) {
      error = true;
    }
  }

  if (args.fasta_sel)
    names.fasta_sel = true;
  else
    names.fasta_sel = false;

  found1 = false;
  if (not error && args.fastq && args.lista) {
    selectfastq(names);
    found1 = true;
    error = true;
  } 
  if (not error && args.fasta && args.lista) {
    selectfasta(names);
    found1 = true;
    error = true;
  } 
  if (not error && args.fastq && args.fastq2fasta) {
    fastq2fasta(names);
    found1 = true;
    error = true;
  } 
  if (not error && args.fastq && args.randn) {
    getrandom(names);
    found1 = true;
    error = true;
  }
 
  if (error && not found1)
    printerror(argv[0]);
  
  return 1;
}

void getrandom(Names& names) 
{
  string line, salida;
  ifstream in;
  unordered_map <int, int> lista;
  unordered_map <int, int>::const_iterator got;
  srand ( unsigned ( time(0) ) );

  
  int total_lines = wc(names.fastq.c_str())/4;

  int count = 0;
  while (count < names.randn) {
    int newn = rand()%total_lines;
    got = lista.find(newn); 
    if ( got == lista.end() ) {
      lista[newn] = 1;
      count++;
    }
  }

  in.open(names.fastq.c_str());

  unsigned char tipo = 0;
  bool mostrar = false;
  while (getline(in,line).good()) {
    if (tipo == 0) {
      mostrar = false;
      got = lista.find(count);
      if ( got != lista.end() ) {
        mostrar = true;
        cout << line << endl;
      }
      count++;
    } else {
      if (mostrar)
        cout << line << endl;
    }
    tipo = (tipo + 1) % 4;
    
  }
  in.close();
}

void selectfasta(Names& names) 
{
  string line, salida;
  ifstream in;
  unordered_map <string, int> lista;
  unordered_map<string, int>::const_iterator got;
  string word;

  in.open(names.lista.c_str());
  bool andpersant = false;
  bool good = getline(in,line).good();
  if (good) {
    if (line.at(0) == '>') {
      istringstream iss(line.substr(1,1000));
      iss >> word;
      lista[word] = 1;
      while (getline(in,line).good()) {
        istringstream iss2(line.substr(1,1000));
        iss2 >> word;
        lista[word] = 1;
      }
    } else {
      istringstream iss(line);
      iss >> word;
      lista[word] = 1;
      while (getline(in,line).good()) {
        istringstream iss2(line);
        iss2 >> word;
        lista[word] = 1;
      }
    }
  }
  in.close();

  in.open(names.fasta.c_str());
  bool mostrar = false;
  while (getline(in,line).good()) {

    if (line.at(0) == '>') {
      mostrar = false;
      istringstream iss(line.substr(1,1000));
      iss >> word;
      got = lista.find(word);
      if ( got != lista.end() ) {
        if (names.fasta_sel) {
          mostrar = true;
          cout << line << endl;
        }
      } else {
        if (not names.fasta_sel) {
          mostrar = true;
	  cout << line << endl;

        }
      }

    } else {
      if (mostrar)
        cout << line << endl;
    }
  }
  in.close();
}

void selectfastq(Names& names) 
{
  
  string line, salida;
  ifstream in;
  unordered_map <string, int> lista;
  unordered_map<string, int>::const_iterator got;

  in.open(names.lista.c_str());
  bool andpersant = false;
  bool good = getline(in,line).good();
  if (good) {
    if (line.at(0) == '@') {
      lista[line.substr(1,1000)] = 1;
      while (getline(in,line).good()) 
        lista[line.substr(1,1000)] = 1;
    } else {
      lista[line] = 1;
      while (getline(in,line).good()) 
        lista[line] = 1;
    }
  }
  in.close();

  bool mostrar_enable = false;
  if (not names.fasta_sel)
    mostrar_enable = true;
  
  in.open(names.fastq.c_str());
  unsigned char tipo = 0;
  bool mostrar = mostrar_enable;
  int title=0;
  while (getline(in,line).good()) {
    if (tipo == 0) {
      mostrar = mostrar_enable;
      istringstream iss(line.substr(1,1000));
      string word;
      iss >> word;
      got = lista.find(word);
      if ( got != lista.end() ) {
        mostrar = not (mostrar_enable);
        //cout << line << endl;
      }

    } 
    //else {
      if (mostrar)
        cout << line << endl;
    //}
    tipo = (tipo + 1) % 4;
  }
  in.close();
}

void fastq2fasta(Names& names) 
{
  string line, salida;
  ifstream in;

  in.open(names.fastq.c_str());
  unsigned char tipo = 0;
  bool mostrar = false;
  while (getline(in,line).good()) {
    if (tipo == 0) 
      cout << ">" << line.substr(1,1000) << endl; 
    if (tipo == 1)
      cout << line << endl;
    
    tipo = (tipo + 1) % 4;
  }
  in.close();
}

int mstrlen(const char arg[]) {
  int out = 0;
  while(arg[out] != 0)
    out++;
  return out;
}

bool IsParam(char arg[],const char comp[])
{
  int len_comp = mstrlen(comp);
  bool equal = true;
  if (mstrlen(arg) == len_comp) {
    for (int i = 0; i < len_comp; i++) {
      if (arg[i] != comp[i]) {
        equal = false;
      }
    }
  } else {
    equal = false;
  }
  return equal;
}

static uintmax_t wc(char const *fname)
{
    static const auto BUFFER_SIZE = 16*1024;
    int fd = open(fname, O_RDONLY);
    if(fd == -1) {
        cerr << "could not open file" << fname << " for reading" << endl;
      return 0;
    }
    /* Advise the kernel of our access pattern.  */
    posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL

    char buf[BUFFER_SIZE + 1];
    uintmax_t lines = 0;

    while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
    {
        if(bytes_read == (size_t)-1){
            cerr << "failed to read " << fname << endl;
            return 0;
        }
        if (!bytes_read)
            break;

        for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
            ++lines;
    }

    return lines;
}


void printerror(const char arg[])
{
    cout << "Select fastQ or fasta reads from list, -list or -random is required" << endl;
    cout << endl;
    cout << arg << endl;
    cout << "  -h                  (Help)" << endl;
    cout << "  -fastq        FILE  (fastq file to select reads from) " << endl;
    cout << "  -list         FILE  (list of reads, fastq or fasta) " << endl;
    cout << "  -random       VAL   (number of random reads to be selected from fastq file) " << endl;
    cout << "  -fastq2fasta        (convert fastq file to fasta)" << endl;
    cout << "  -fasta        FILE  (fasta file to select reads from) " << endl;
    cout << "  -fasta_sel          (from fasta file select reads in -list, if not flag, reads not in list are selected)" << endl;
    cout << ""<< endl;

}
