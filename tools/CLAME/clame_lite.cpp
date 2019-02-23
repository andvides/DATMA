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

#include "clame_lite.h"


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
        if ( IsParam(argv[argum],"-e") ) //edges constrains
        {
            argum ++;
            if (argum < argc)  
            {
                args->edges = true;
                std::string::size_type sz;     
                parameters->edges = std::stof (argv[argum],&sz);
                //cout<<"edges:"<<parameters->edges<<endl;
           } 
            else 
                error = true;
         
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-tol") ) //MAD tol
        {
            argum ++;
            if (argum < argc)  
            {
                args->tol = true;
                std::string::size_type sz;     
                parameters->tol = std::stof (argv[argum],&sz);
                //cout<<"edges:"<<parameters->edges<<endl;
           } 
            else 
                error = true;
         
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-w") ) // windows offset to start the alignemnt
        {
            argum ++;
            args->w = true;
            found1 = true;
            if (argum < argc)  
                parameters->w = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-ld") ) // minimal number of total edges
        {
            argum ++;
            args->ld = true;
            found1 = true;
            if (argum < argc)  
                parameters->ld = atoi(argv[argum]);
            argum ++;
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
        if ( IsParam(argv[argum],"-print") ) //print to file
        {
            argum ++;
            parameters->enablePrint = true;
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
         if ( IsParam(argv[argum],"-fastq") ) 
        {
            argum ++;
            args->fastq = true;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-fr") ) //to print the help
        {
            argum ++;
            if (argum < argc)  
            {
                args->forward_reverse = true;
                parameters->forward_reverse = argv[argum];
            } 
            else 
                error = true;
                       
            argum ++;
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

bool alignemnt(Names *names,Parameters *parameters,vector<string> *bases,string *fasta, vector<uint32_t> *index, int* queryList, vector<int>* MatrixList)
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
           
            for (auto it=locations.begin(); it<locations.end(); it++)
                index->push_back(*it);
        }
    }
    else//buidl suffix tree
    {
    
        index_file   = names->outputFile+".fm9";
        construct_im(fm_index,*fasta,1); 
        *fasta="";
        if(parameters->enablePrint)
            store_to_file(fm_index, index_file); // save it
    }
    
    //search
    int cpus=parameters->numThreads;
    int numberOFreads=bases->size();
    bool runningError=false;
    
    string typeSearch=parameters->forward_reverse;
    
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
                int totalSearch=4;
                if (query.length()>=w+parameters->query_size)
                {
                    string reverseQuery=reverse(query);
                    string queryF1=query.substr(w,parameters->query_size);//forward init
                    string queryF2=query.substr((query.length()-w-parameters->query_size),parameters->query_size);;//forward end
                    string queryR1=reverseQuery.substr(w,parameters->query_size);//reverse init
                    string queryR2=reverseQuery.substr((reverseQuery.length()-w-parameters->query_size),parameters->query_size);//reverse end
                    string Queries[4]={queryF1,queryF2,queryR1,queryR2};
                    string orientations[4]={"5a'","5b'","3a'","3b'"};
                    
                    if(typeSearch=="f")
                    {
                        totalSearch=2;
                        Queries[0]=queryF1;
                        Queries[1]=queryF2;
                    
                    }
                    else if (typeSearch=="r")
                    {
                        totalSearch=2;
                        Queries[0]=queryR1;
                        Queries[1]=queryR2;
                    }
                    else if (typeSearch=="f1r1")
                    {
                        totalSearch=2;
                        Queries[0]=queryF1;
                        Queries[1]=queryR1;
                    }
                    else if (typeSearch=="f1r2")
                    {
                        totalSearch=2;
                        Queries[0]=queryF1;
                        Queries[1]=queryR2;
                    }
                    else if (typeSearch=="f2r1")
                    {
                        totalSearch=2;
                        Queries[0]=queryF2;
                        Queries[1]=queryR1;
                    }
                    else if (typeSearch=="f2r2")
                    {
                        totalSearch=2;
                        Queries[0]=queryF2;
                        Queries[1]=queryR2;
                    }



                    queryList[i]=0;// Generate the nodes of the tree and init as unmarked read
                    for (int q=0;q<totalSearch;q++)
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
                                MatrixList[subjectID].push_back(i);
                            }
                        }
                    }
                }
                else
                    queryList[i]=10000;// Read too short
               
           
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
        if(parameters->enablePrint)
            printResult(names,numberOFreads,MatrixList);

        return runningError;
    }
}

void printResult(Names *names, int numberOFreads, vector<int> *MatrixList)
{
    ofstream myfile, myfile2;
    string nameFile=names->outputFile+".result";//ss.str();
    string nameFile2=names->outputFile+".links";//ss.str();
    myfile.open(nameFile.c_str());
    myfile2.open(nameFile2.c_str());

    for(int k=0;k<numberOFreads;k++)
    {
        //print links
        myfile2<<k<<"\t"<<MatrixList[k].size()<<endl;
        
        //print result
        myfile<<k<<"\t";
        for (auto it=MatrixList[k].begin(); it<MatrixList[k].end(); it++)
            myfile<<*it<<"\t";
        myfile<<"\n";
    }
    myfile.close();
    myfile2.close();
}


void binningPrint(Names *names, Parameters *parameters, vector<string> *title, int *queryList, vector<int> *MatrixList, int numberOFreads, vector<string> *bases,vector<string> *qual)
{
    
    /*std::vector<int> edges_vect;
    std::stringstream ss(parameters->edges);
    int i;

    while (ss >> i)
    {
        edges_vect.push_back(i);

        if (ss.peek() == ',')
            ss.ignore();
    }
    std::sort (edges_vect.begin(), edges_vect.end()); */
    
    
    int cpus=parameters->numThreads;
    std::vector<string>::iterator ptrBases= bases->begin();
    std::vector<string>::iterator ptrTitle= title->begin();
    std::vector<string>::iterator ptrQual=  qual->begin();
    
    
    
    //1.Number of links by read
    #pragma omp parallel num_threads(cpus)
    {
        #pragma omp for schedule(runtime)
        for(int k=0; k<numberOFreads;k++) 
            if(MatrixList[k].size()<=parameters->ld )
                queryList[k]=1000; //marked as  visited 
    }
    
    

    
    ofstream myfile;
    string nameFile=names->outputFile+".binning";
    myfile.open(nameFile.c_str());
    //myfile<<"Bin\tld\tlu\tsize\tmean\tstdv\tmean1/3\tcv1\tmedian\tmode\tmad"<<endl;
    myfile<<"@bin\tsize\tmean\tstd\tmedian\tMAD\tp=3std/mean\tdist2mean\tOutliers"<<endl;


    int numBin=0; 
    
    //Genera los bins
    string ext=".fasta";
    if(parameters->fastq)
        ext=".fastq";
        
      
    int *stack = new int[numberOFreads] ();
    int *put=stack;
    int *get=stack;

    for(int i=0; i<numberOFreads;i++)
    {
        if(queryList[i]==0) //not visited
        {
            queryList[i]=1;
            *(put++)=i;
            while(get!=put)
            {
                int q=*get++;
                for( vector<int>::iterator j=MatrixList[q].begin(); j!=MatrixList[q].end(); ++j)
                {
                    if(queryList[*j]==0) //not visited
                    {
                        queryList[*j]=1;
                        *(put++)=*j;
                    }
                }
            }

            //BIN static_parameters
            int size=get-stack;
            int maximun, minimum;
            float cut=parameters->edges;
            float *stadistic = new float[7] (); 
            if(size>=parameters->sizeBin)
            {
                
                static_parameters(stack, size, MatrixList, queryList, numBin, &maximun, &minimum,cut, stadistic);
                myfile<<"@bin"<<numBin<<"\t"<<size<<"\t"<<stadistic[0]<<"\t"<<stadistic[1]<<"\t"<<stadistic[3]<<"\t"<<stadistic[5]<<"\t"<<stadistic[6]<<"\t"<<cut<<"\t"<<minimum<<"\t"<<maximun<<endl;
                //cout<<"@bin"<<numBin<<"\t"<<size<<"\trange:\t"<<minimum<<"\t"<<maximun<<"\t"<<stadistic[0]<<"\t"<<stadistic[1]<<"\t"<<stadistic[6]<<"\t"<<stadistic[2]<<endl;

                //numBin++;
                int *stack2 = new int[size] ();
                int itertation=0;
                //float error=abs(stadistic[6]-1);
                float tol=parameters->tol;//2.0;////0.5;
                float red=0.25;
                while((abs(stadistic[6]-1.0)>tol) && (cut>0 ) && (size>1) && stadistic[6]>1.0)//while(stadistic[6]>1 && cut>0)
                {
        
                    //Number of links by read
                    //#pragma omp parallel num_threads(cpus)
                    //{
                        //#pragma omp for schedule(runtime)
                        for(int k=0; k<size;k++) 
                        {
                            int s=*(stack+k);
                            if(MatrixList[s].size()>maximun || MatrixList[s].size()<minimum)
                                queryList[s]=0; //marked as  non taked queryList[s]=3
                            else
                                queryList[s]=1;
                        }
                        //cout<<"deleted "<<p2<<endl;
                    //}
                
                    //Redefine the binning
                    put=stack2;
                    get=stack2;
                    for(int k=0; k<size;k++)
                    {
                        int s=*(stack+k);
                        if(queryList[s]==1) //no binned
                        {
                            queryList[s]=2;
                            *(put++)=s;
                            while(get!=put)
                            {
                                int q=*get++;
                                for( vector<int>::iterator j=MatrixList[q].begin(); j!=MatrixList[q].end(); ++j)
                                {
                                    if(queryList[*j]==1) //not binned
                                    {
                                        queryList[*j]=2;
                                        *(put++)=*j;
                                    }
                                }
                            }
                        }
                    }
                
                    size=get-stack2;
                    static_parameters(stack2, size, MatrixList, queryList, numBin, &maximun, &minimum,cut, stadistic);
                    //myfile<<">bin"<<numBin<<"\t"<<size<<"\trange:\t"<<minimum<<"\t"<<maximun<<endl;
                    //p=3*stadistic[1]/stadistic[0];
                    myfile<<itertation<<"_bin"<<numBin<<"\t"<<size<<"\t"<<stadistic[0]<<"\t"<<stadistic[1]<<"\t"<<stadistic[3]<<"\t"<<stadistic[5]<<"\t"<<stadistic[6]<<"\t"<<cut<<"\t"<<minimum<<"\t"<<maximun<<endl;
                    
                    itertation++;
                    cut-=red;
                }
                
                if(cut==parameters->edges) //no statistic control necessary
                    stack2=stack;
                
                //Print the Bin
                ofstream myfile2;
                stringstream ss2;
                ss2 << numBin;
                string strNumBin = ss2.str();
                string nameFile2=names->outputFile+"_"+strNumBin+ext;
                myfile2.open(nameFile2.c_str());

                while(get!=stack2)
                {
                    --get;
                    //myfile<<*get<<"\t"<<numBin<<endl;
                    string base=*(ptrBases+*(get));
                    string name=*(ptrTitle+*(get));
                    myfile2<<name<<endl;
                    myfile2<<base<<endl;
                    if(parameters->fastq)
                        {string quality=*(ptrQual+*(get)); myfile2<<'+'<<endl, myfile2<<quality<<endl;}
                
                }
                                
                myfile2.close();
                numBin++;
                
            }
            
            put=stack; //pointers restart
            get=stack; 
        }
    }
    myfile.close();

}
void static_parameters(int *stack, int size, vector<int> *MatrixList, int *queryList, int numBin, int *maximun, int *minimum, float cut, float *stadistic) 
{
    float mean=0.0;
    float std,std1=0.0;
    float cv;      
    float median;
    float mode;
    
    vector <float> v,v2;
    vector <int> v3;
    
    //mean
    float links=0.0;
    for( int p=0;p<size;p++)
    {
        int edges=MatrixList[*(stack+p)].size();
        links+=edges;
        v.push_back(edges);
    }
    mean=(float)(links/size);
    stadistic[0]=mean;
    
    //std
    for( int p=0;p<size;p++)
    {
        int edges=MatrixList[*(stack+p)].size();
        std1+=pow(edges - mean, 2);
    }
    std=sqrt(std1/size);
    stadistic[1]=std;
    
    //variation coefficient
    cv=100*std/mean;          
    stadistic[2]=cv;
    
    //median
    sort(v.begin(), v.end());
    median = size % 2 ? v[size/2] : (v[size/2-1] + v[size / 2] / 2);
    stadistic[3]=median;
    
    //mode
    float number = v[0];
    mode = number;
    int count = 0;
    int countMode = 1;

    for (int i=0; i<size; i++)
    {
        if (v[i] == number) 
            count++;
        else
        {
            if (count > countMode) 
            {
                  countMode = count;
                  mode = number;
            }
           count = 1;
           number = v[i];
        }
    }
    stadistic[4]=mode;
    
    //mad 
    for (int i=0; i<size; i++)
    {
        float aux=pow(v[i]-median, 2);
        v2.push_back(sqrt(aux));
    }
    sort(v2.begin(), v2.end());
    float median2 = size % 2 ? v2[size/2] : (v2[size/2-1] + v2[size / 2] / 2);
    float mad = 1.4826*median2;
    stadistic[5]=mad;
    
    //cutoff
    for (int i=0; i<size; i++)
    {
        int edges=MatrixList[*(stack+i)].size();
        float aux=sqrt(pow(edges-median, 2));  
        float cutoff= aux/mad;
        //cout<<">bin"<<numBin<<"\t"<<size<<"\t"<<*(stack+i)<<"\t"<<edges<<"\t"<<cutoff<<endl;
        if(cutoff<cut)
             v3.push_back(edges);
    }
    
    //sort(v3.begin(), v3.end());
    int max=v3[0]; int min=v3[0];
    for (std::vector<int>::iterator it = v3.begin() ; it != v3.end(); ++it)
    {
        if(*it>max)
            max=*it;
        
        if (*it<min)
            min=*it;
            
    }
    *(maximun)=max;
    *(minimum)=min;
    stadistic[6]=3*stadistic[1]/stadistic[0];

    //cout<<"@bin"<<numBin<<"\tsize:\t"<<v3.size()<<"\trange:\t"<<*(v3.begin())<<"\t"<<*(v3.end())<<"\t"<<mini<<"\t"<<max<<endl;
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
    cout << "CLAME:'Clasificador Metagenomico'" << endl;
    cout << "version 2.3 Febrary 2017"<<endl;
    cout << "Authors"<<endl;
    cout << "Benavides A, Alzate JF and Cabarcas F"<<endl;
    cout << endl;
    cout << arg << endl;
    cout << "  -h\t\t\t(Help)" << endl;
    cout << "  -b minimum number of bases to take an alignment (default 70) " << endl;
    cout << "  -e normal cut point for notmal test (default 3) " << endl;
    cout << "  -tol MAD error tolerance (default 0.5) " << endl;
    cout << "  -fm9 Load fm9 file  " << endl;
    cout << "  -fastq input file is in a fastq format  " << endl;
    cout << "  -fr disable forward or reverse search in the alignemnt (f=forwar, r=reverse, fr=forwar and reverse (default)) " << endl;
    cout << "  -multiFasta\t\tFILE  with all the reads " << endl;
    cout << "  -ld minimal nnumber of total edges  (default 2) " << endl;
    cout << "  -nt number of threads to use (default 1) " << endl;
    cout << "  -output name for the output-file  if print option was selected (default output)" << endl;
    cout << "  -print enable print output to file (default false) " << endl;
    cout << "  -sizeBin minimum number of reads to report a bin (default 1000) " << endl;
    cout << "  -w windows offset to start the alignemnt (default 0) " << endl;
    cout << ""<< endl;

}
