/*---------------------------------------------------------------
 *
 *   CLAME:"mapping"
 *	
 *   
 *  
 ---------------------------------------------------------------*/

#include "mapping_lite.h"

bool readArguments(int argc,char *argv[], Args *args,Names *names, Parameters *parameters)
{
    bool   error = false, found1=false, mandatory=false;
    int    argum = 1;
    while (argum < argc && not error) 
    {
        found1 = false; 
        if ( IsParam(argv[argum],"-multiFasta") ) //multifasta file
        {
            argum ++;
            if (argum < argc)  
            {
                args->multiFasta = true;
                names->multiFasta = argv[argum];
                mandatory=true;
            } 
            else 
                error = true;
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-fm9") ) //FM) reference
        {
            argum ++;
            if (argum < argc)  
            {
                args->fm9 = true;
                names->fm9 = argv[argum];
            } 
            else 
                error = true;
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-lu") ) //superior cut
        {
            argum ++;
            if (argum < argc)  
            {
                args->lu = true;
            parameters->lu = atoi(argv[argum]);
            } 
            else 
                error = true;
         
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-ld") ) //inferior cut
        {
            argum ++;
            if (argum < argc)  
            {
                args->ld = true;
                parameters->ld = atoi(argv[argum]);
            } 
            else 
                error = true;
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-output") ) //output file
        {
            argum ++;
            if (argum < argc)  
            {
                args->outputFile = true;
                names->outputFile = argv[argum];
            } 
            else 
                error = true;
        
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-nt") ) // number of cpus
        {
            argum ++;
            args->numT = true;
            found1 = true;
            if (argum < argc)  
                parameters->numThreads = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-b") ) //minimal base lenght
        {
            argum ++;
            args->bases_Threshold = true;
            found1 = true;
            if (argum < argc)  
                parameters->query_size = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-w") ) //minimal reads set to report a bin
        {
            argum ++;
            if (argum < argc)  
            {
                args->w = true;
                parameters->w = atoi(argv[argum]);
            } 
            else 
                error = true;
                       
            argum ++;
            found1 = true;
            continue;
        }       
        if ( IsParam(argv[argum],"-sizeBin") ) //minimal reads set to report a bin
        {
            argum ++;
            if (argum < argc)  
            {
                args->sizeBin = true;
                parameters->sizeBin = atoi(argv[argum]);
            } 
            else 
                error = true;
                       
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-print") ) 
        {
            argum ++;
            args->print = true;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-fastq") ) 
        {
            argum ++;
            args->fastq = true;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-h") ) //to print the help
        {
            error = true;
            found1 = false; 
            continue;
        }
        if (not found1)
            error = true;
    }//end of the while

    if (error && not found1 ||  argc < 2  || not mandatory)
        return false;
    else
        return true;
}
bool readFasta(Names *names, vector<string> *bases,string *fasta,vector<uint32_t> *index, vector<string> *title)
{
    
    //to print index
    ofstream myfile;
    string nameFile=names->outputFile+".index";
    myfile.open(nameFile.c_str());

    ifstream infile;
    infile.open(names->multiFasta.c_str());
    
    if (infile.is_open()) 
    {
        string word, name="", read="";
        bool firstTime=true;
        int indexBases=0;
        int indexTitle=0;
        uint32_t n=0;
        index->push_back(indexBases);

        while (getline(infile,word))
        {
            if (word.find(">") != string::npos) 
            {
                if (firstTime)
                    firstTime=false;
                else
                {
                    myfile<<name<<"\t"<<indexTitle++<<endl; 
                    title->push_back (name);
                    bases->push_back(read);
                    read+='$';
                    n = read.length();
                    *fasta+=read;
                    indexBases+=n;
                    index->push_back(indexBases);
                    read="";
                }
                int found = word.find(' ');
                name=word.substr (0,found);
            }
            else
                read=read+word;
        }
        myfile<<name<<"\t"<<indexTitle++<<endl; 
        title->push_back (name);
        bases->push_back(read);
        read+='$';
        n = read.length();
        *fasta+=read;
        indexBases+=n;
        index->push_back(indexBases);
        read="";
        infile.close();
        
        myfile.close();
        return false;
    }
    else
    {
        #ifdef debug
            std::cout << "stage1: Error opening fasta file";
        #endif
        return true;
    }
}

bool readFastQFile(Names *names,  vector<string>* bases, string* fasta, vector<uint32_t> *index, vector<string> *title,  vector<string> *qual)
{
    //to print the index mandatory file
    ofstream myfile;
    string nameFile=names->outputFile+".index";//ss.str();
    myfile.open(nameFile.c_str());

    int indexTitle=0;
    int indexBases=0;
    uint32_t n=0;

    //name file string word
    ifstream infile;
    infile.open(names->multiFasta.c_str());
    typedef enum {s0, s1,s2} STATES; //to read multifasta files
    STATES state=s0;

    if (infile.is_open()) 
    {
        string word, name="", read="", quality="";
        bool firstTime=true;
        int countLine=0;
        index->push_back(indexBases);
        while (getline(infile,word))
        {
            switch ( state ) 
            {
                case s0:
                    if (word[0]=='@') 
                    {    
                        int found=word.find(' '); 
                        name=word.substr(0,found);
                        countLine=0; state=s1; myfile<<name<<"\t"<<indexTitle++<<endl; //print Read0 0...;
                        title->push_back (name);
                    }
                    break;
                case s1:
                    if (word[0]=='+') 
                        state=s2;
                    else
                        read=read+word, countLine++;
                    break;
                case s2:
                    quality=quality+word;
                    if(--countLine==0)
                    {
                        qual->push_back(quality);
                        bases->push_back(read);
                        read+='$';
                        n = read.length();
                        indexBases+=n;
                        index->push_back(indexBases);
                        *fasta+=read;
                        read="";
                        quality="";
                        state=s0;
                    }
                    break;
            }
        }
        infile.close();
        myfile.close();
        return false;
    }
    else
    {
        #ifdef debug
            std::cout << "stage1: Error opening fastq file";
        #endif
        return true;
    }
}


bool alignemnt(Names *names,Parameters *parameters,vector<string> *bases,string *fasta, vector<uint32_t> *index, int* queryList, vector<int>* MatrixList, vector<string> *title)
{
       
    std::vector<string>::iterator ptrBases= bases->begin();
   
    //FM-Index
    csa_wt<> fm_index;
    string index_file =""; 
    if(parameters->loadFM9)//load from file
    {

        index_file   = names->fm9;
        if (!load_from_file(fm_index, index_file)) 
        {
            #ifdef debug
                cout << "No index "<<index_file<< "error" << endl;
            #endif
            return true;
        }
        else
        {
            index->clear();
            index->push_back(0);
            string querySearch="$";
            auto locations = locate(fm_index, querySearch.begin(), querySearch.begin()+1);
            sort(locations.begin(), locations.end());
            int i=0;
            for (auto it=locations.begin(); it<locations.end(); it++)
	    {
                index->push_back(*it);
                #ifdef debug
		    cout<<i++<<"\t"<<*it<<endl;
            	#endif
	    }
        }
    }
    else//buidl suffix tree
    {
    
            #ifdef debug
                cout << "No index file found" << endl;
            #endif
            return true;
    }
    
    //search
    int cpus=parameters->numThreads;
    int numberOFreads=bases->size();
    bool runningError=false;
    
    if(runningError==false)
    {
        #pragma omp parallel num_threads(cpus)
        {
            #pragma omp for schedule(runtime) 
            for(int i=0; i<numberOFreads;i++)//process all parts except the last one because it can get a diffent size
            {
                uint32_t subjectID=0;
                string query=*(ptrBases+i);
                int w=parameters->w;
                if (query.length()>=w+parameters->query_size)
                {
                    string reverseQuery=reverse(query);
                    
                    string queryF1=query.substr(w,parameters->query_size);//forward init
                    string queryF2=query.substr((query.length()-w-parameters->query_size),parameters->query_size);;//forward end
                    string queryR1=reverseQuery.substr(w,parameters->query_size);//reverse init
                    string queryR2=reverseQuery.substr((reverseQuery.length()-w-parameters->query_size),parameters->query_size);//reverse end
    
                    string Queries[4]={queryF1,queryF2,queryR1,queryR2};
                    string orientations[4]={"5a'","5b'","3a'","3b'"};
	
                    queryList[i]=0;// Generate the nodes of the tree and init as unmarked read
                    for (int q=0;q<4;q++)
                    {
                        string querySearch=Queries[q];
                        uint32_t sizeq = querySearch.length();
                        auto locations = locate(fm_index, querySearch.begin(), querySearch.begin()+sizeq);
                        for (auto it=locations.begin(); it<locations.end(); it++)
                        {
                            subjectID=indx2Loc(*it,index);
                            #pragma omp critical
                            {
                                MatrixList[i].push_back(subjectID);
                                //MatrixList[subjectID].push_back(i);
                            }
                        }
                    }
                }
            }
        }//end pragma
        
        //sort and deleted duplicated items to get the number of links
        for(int k=0;k<numberOFreads;k++)
        {
            std::sort (MatrixList[k].begin(), MatrixList[k].end());           
            std::vector<int>::iterator it=std::unique (MatrixList[k].begin(), MatrixList[k].end()); 
            MatrixList[k].resize( std::distance(MatrixList[k].begin(),it) );
        }
        
        //print the result
        printResult(names,parameters,numberOFreads,MatrixList,title);
            
        return runningError;
    }
}

void printResult(Names *names, Parameters *parameters,int numberOFreads, vector<int> *MatrixList, vector<string> *title)
{
    std::vector<string>::iterator ptrTitle= title->begin();

    ofstream myfile, myfile2;
    string nameFile2=names->outputFile+".links";//ss.str();
    myfile2.open(nameFile2.c_str());

    if (parameters->enablePrint)
    {
	//ofstream myfile;
    	string nameFile=names->outputFile+".result";//ss.str();
    	myfile.open(nameFile.c_str());
    }

    for(int k=0;k<numberOFreads;k++)
    {
        if(MatrixList[k].size()>=parameters->ld && MatrixList[k].size()<=parameters->lu )
	{
		//print links
            	myfile2<<*(ptrTitle+k)<<"\t"<<MatrixList[k].size()<<endl;
        
        	//print result
		if (parameters->enablePrint)
		{
	 		myfile<<k<<"\t";
        		for (auto it=MatrixList[k].begin(); it<MatrixList[k].end(); it++)
            			myfile<<*it<<"\t";
        		myfile<<"\n";
		}
	}
    }
    myfile2.close();
    if (parameters->enablePrint)
        myfile.close();
}

void binningPrint(Names *names, Parameters *parameters, vector<string> *title, int *queryList, vector<int> *MatrixList, int numberOFreads, vector<string> *bases)
{
    
}

uint32_t indx2Loc(uint32_t locs, vector<uint32_t> *index)
{
        std::vector<uint32_t>::iterator up;
        up= std::upper_bound (index->begin(), index->end(), locs); 
        return (up-1-index->begin());
}


string reverse(string str)
{
    // This table is used to transform nucleotide letters into reverse. //A=T, C=G, G=C, T=A, Other=N
    static const char  nt_table_reverse[128] = {
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'A', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'A', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
        };
  
        string rs="";
        for (std::string::reverse_iterator rit=str.rbegin(); rit!=str.rend(); ++rit)
        {
            char base=nt_table_reverse[*rit]; 
            rs=rs+base;
        }
        //cout<<str<<"\t"<<rs<<endl;
        return rs;


}

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

void printerror(const char arg[])
{
    cout << "CLAME: mapping to reference" << endl;
    cout << endl;
    cout << arg << endl;
    cout << "  -h\t\t\t(Help)" << endl;
    cout << "  -b minimum number of bases to take an alignment (default 20) " << endl;
    cout << "  -fm9 Load fm9 file  " << endl;
    cout << "  -fastq input file is in a fastq format  " << endl;
    cout << "  -ld minimun number of links (default 1) " << endl;
    cout << "  -lu maximun number of links (default 10000) " << endl;
    cout << "  -multiFasta\t\tFILE  with all the reads " << endl;
    cout << "  -nt number of threads to use (default 1) " << endl;
    cout << "  -output name for the output-file  if print option was selected (default output)" << endl;
    cout << "  -print print the result file (default false)" <<endl;
    cout << "  -w windows offset to start the alignemnt (default 0) "<<endl;
    cout << ""<< endl;

}




