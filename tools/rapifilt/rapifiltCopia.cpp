/* L I C E N S E *************************************************************/

/*
    Copyright (C) 2009, 2010 Indraniel Das <indraniel@gmail.com>
                             and Washington University in St. Louis

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, see <http://www.gnu.org/licenses/>
*/

/* I N C L U D E S ***********************************************************/
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
#include "sff.h"


/* D E F I N E S *************************************************************/
#define VERSION "V0"
#define PRG_NAME "rapifilt"
#define FASTQ_FILENAME_MAX_LENGTH 1024
#define input_FILENAME_MAX_LENGTH 1024

using namespace std;


/* P R O T O T Y P E S *******************************************************/
void process_options(int argc, char *argv[], struct Args * args, struct Names * names);
int mstrlen(const char arg[]);
bool IsParam(char arg[],const char comp[]);
void help_message(void);
void process_illumina_to_fastq(struct Args * args, struct Names * names);
void process_sff_to_fastq(struct Args * args, struct Names * names);
void construct_fastq_trimmed(FILE *fp, char *name, char *bases, uint8_t *quality, int nbases, struct Args * args, struct Names * names);

/* G L O B A L S *************************************************************/
struct Args 
{
	bool outputFile;	/* Output file name */ 
	bool illuminaFile; 	/* Input illumina file name */ 
	bool sffFile; 		/* Input illumina file name */ 
	bool enableFasta; 	/* selec between fasta(1) and fastq output(0)*/
	int left_cut; 		/* left cut value */
	int right_cut; 		/* right cut value */
	int windows_size;	/* windows size to check */
	int min_len;		/* Filter sequence shorter than min_len*/	
	int max_len;		/* Filter sequence larger than max_len*/	
};
struct Names {string outputFile; string illuminaFile1; string illuminaFile2; string sffFile;};

/* M A I N *******************************************************************/
int main(int argc, char *argv[]) 
{
	Args   args = {false, false, false, false, 0, 0, 1, 0,5000};
  	Names  names;

  	process_options(argc, argv, &args, &names);

	if(args.illuminaFile)
		process_illumina_to_fastq(&args, &names);
    	else if	(args.sffFile)	
    		process_sff_to_fastq(&args, &names);
	else
		help_message(),	exit(0);
	return 0;
}

/* F U N C T I O N S *********************************************************/
void process_illumina_to_fastq(struct Args * args, struct Names * names)
{
	//output file
	FILE *Outfastq_fp1, *Outfastq_fp2, *Singleton_fp1, *Singleton_fp2, *bad_fp1, *bad_fp2, *fp1, *fp2;
	string ext=".fastq"; 
	char firstC='@';
	if(args->enableFasta)
		firstC='>',ext=".fasta";

	if ( args->outputFile ) 
    	{
		string name5=names->outputFile+"_1"+ext;
		string name3=names->outputFile+"_2"+ext;
		string single1=names->outputFile+"_1Singleton"+ext;
		string single2=names->outputFile+"_2Singleton"+ext;
		string bad1=names->outputFile+"_1bad"+ext;
		string bad2=names->outputFile+"_2bad"+ext;

		char myArray1[name5.size()+1];
		char myArray2[name3.size()+1];		
		char myArray3[single1.size()+1];
		char myArray4[single2.size()+1];
		char myArray5[bad1.size()+1];
		char myArray6[bad2.size()+1];
    		strcpy(myArray1, name5.c_str());
    		strcpy(myArray2, name3.c_str());
		strcpy(myArray3, single1.c_str());
    		strcpy(myArray4, single2.c_str());
		strcpy(myArray5, bad1.c_str());
		strcpy(myArray6, bad2.c_str());
        	if ( (Outfastq_fp1 	= fopen(myArray1, "w")) == NULL || (Outfastq_fp2 	= fopen(myArray2, "w"))	== NULL ||
		     (Singleton_fp1 	= fopen(myArray3, "w")) == NULL || (Singleton_fp2 	= fopen(myArray4, "w"))	== NULL ||
		     (bad_fp1 		= fopen(myArray5, "w")) == NULL || (bad_fp2 		= fopen(myArray6, "w"))	== NULL) 	
        		cout<<"Could not open files for writing "<<names->outputFile<<endl, exit(1);
    	}

	//to take each name file
	/*string name [2];
	char * dup = strdup(names->illuminaFile.c_str());
	char * token = strtok(dup, ",");
	int i=0;
	while(token != NULL)
        	name[i++]=(string(token)),token = strtok(NULL, ",");
    	free(dup);
	char * file1 = strdup(name[0].c_str());
	char * file2 = strdup(name[1].c_str());*/

	char * file1 = strdup(names->illuminaFile1.c_str());
	char * file2 = strdup(names->illuminaFile2.c_str());
	//cout<<"Input Files: "<<names->illuminaFile1<<", "<<names->illuminaFile2<<endl;
	
	//read File
	ifstream infile1, infile2;
    	infile1.open(file1);
	infile2.open(file2);

	typedef enum {s0, s1,s2, s3, s4,s5,s6,s7,s8,s9,s10, s11} STATES; //to read multifasta files
        STATES state=s0;
	
	bool next=true, save=true;
	//char word1='\0', word2='\0'; 
	string name1="", name2="", aux1="", aux2=""; 
	vector<char> bases1, bases2;
	vector<char> quality1, quality2;
	vector<unsigned int> quality1_int, quality2_int;

	int sum=0, k=0;
	int bouderLeft1=0, bouderRight1=0, finalSize1=0, size1=0;
	int bouderLeft2=0, bouderRight2=0, finalSize2=0, size2=0;
	std::vector<unsigned int>::iterator it1b, it1e;
	std::vector<unsigned int>::iterator it2b, it2e;

	while (next)
	{
		char word1='\0', word2='\0'; 
		//cout<<state<<endl;
		switch ( state ) 
		{
			case s0:
				next=false;
				while(word1!='\n')
				{
					if(infile1.get(word1)){
						next=true;
						if(int(word1)!=32 && save) 
							name1+=word1;
						else
							aux1+=word1, save=false;}
					else
						break;
				}
				state=s1;
				save=true;
				break;
			case s1:
				next=false;
				while(word2!='\n')
				{	
					if(infile2.get(word2)){
						next=true;
						if(int(word2)!=32 && save) 
							name2+=word2;
						else
							aux2+=word2, save=false;}
					else
						break;
				}
				state=s2;
				save=true;
				break;
			case s2:
				state=s3;
				if(name1!=name2)
					cout<<"Error files are not equals"<<endl, next=false;
				break;
			case s3:
				while(infile1.get(word1) && word1!='\n')
					bases1.push_back(word1);
				state=s4;
				break;
			case s4:
				while(infile2.get(word2) && word2!='\n')
					bases2.push_back(word2);
				state=s5;
				break;
			case s5:
				while(word1!='\n')
					infile1.get(word1);
				while(word2!='\n')
					infile2.get(word2);
				state=s6;
				break;
			case s6:
				while(infile1.get(word1) && word1!='\n')
					quality1.push_back(word1),quality1_int.push_back((unsigned int)word1-33) ;
				state=s7;
				break;
			case s7:
				while(infile2.get(word2) && word2!='\n')
					quality2.push_back(word2),quality2_int.push_back((unsigned int)word2-33) ;;
				state=s8;
				break;
			case s8:
				//cut by quality
				//left	
				size1=quality1_int.size();
				it1b = quality1_int.begin();
				it1e = quality1_int.end();
				bouderLeft1=0;
    				while (it1b!=it1e) 
    				{
        				sum=0;k=0;
					while((bouderLeft1+k)<size1 && k<args->windows_size && *(it1b+k)>=args->left_cut )
						sum++,k++;

					if(sum==args->windows_size)
						break;
					bouderLeft1++;it1b++;
    				}	
    				//right
				it1b = quality1_int.begin();
				it1e = quality1_int.end();
				bouderRight1=0;
    				while (it1e!=it1b) 
    				{
        				it1e--;sum=0;k=0;
					while((size1-bouderRight1-k)>0 && k<args->windows_size && *(it1e-k)>=args->right_cut)
						sum++,k++;

					if(sum==args->windows_size)
						break;
					bouderRight1++;
    				}

				finalSize1=quality1_int.size()-	bouderLeft1-bouderRight1;
				state=s9;
				break;
			case s9:
				//cut by quality
				//left
				size2=quality2_int.size();
				it2b = quality2_int.begin();
				it2e = quality2_int.end();
				bouderLeft2=0;
    				while (it2b!=it2e) 
    				{
        				sum=0;k=0;
					while((bouderLeft2+k)<size2 && k<args->windows_size && *(it2b+k)>=args->left_cut)
						sum++,k++;

					if(sum==args->windows_size)
						break;
					bouderLeft2++;it2b++;
    				}	
    				//right
				it2b = quality2_int.begin();
				it2e = quality2_int.end();
				bouderRight2=0;
    				while (it2e!=it2b) 
    				{
        				it2e--;sum=0;k=0;
					while((size2-bouderRight2-k)>0 && k<args->windows_size && *(it2e-k)>=args->right_cut)
						sum++,k++;

					if(sum==args->windows_size)
						break;
					bouderRight2++;
    				}

				finalSize2=quality2_int.size()-	bouderLeft2-bouderRight2;
				state=s10;
				break;
			case s10:
				if ( !args->outputFile ) 
	        			fp1 = stdout, fp2 = stdout;
				else
				{
					if((finalSize1>=args->min_len && finalSize1<=args->max_len) && (finalSize2>=args->min_len && finalSize2<=args->max_len)) //11
						fp1=Outfastq_fp1, fp2=Outfastq_fp2;
					else if((finalSize1>=args->min_len && finalSize1<=args->max_len) &&  (finalSize2<args->min_len || finalSize2>args->max_len))//10
						fp1=Singleton_fp1, bouderLeft2=0, bouderRight2=0, fp2=bad_fp2;
					else if((finalSize1<args->min_len || finalSize1>args->max_len) &&  (finalSize2>=args->min_len && finalSize2<=args->max_len))//01
						bouderLeft1=0,bouderRight1=0,fp1=bad_fp1, fp2=Singleton_fp2;
					else //to bad files (next version)
						bouderLeft1=0,bouderRight1=0,fp1=bad_fp1,bouderLeft2=0,bouderRight2=0,fp2=bad_fp2;
				}
				state=s11;
				break;
			case s11:
				name1.erase (0,1);
				//print name
				fprintf(fp1,"%c%s%s",firstC,name1.c_str(),aux1.c_str());//cout<<firstC<<name<<aux1;
				//print bases
				for (std::vector<char>::iterator it = bases1.begin()+bouderLeft1; it != bases1.end()-bouderRight1; ++it)
    					fprintf(fp1, "%c", *it);//std::cout << *it;
				fprintf(fp1, "\n");//cout<<endl;
				//print qualities
				if(args->enableFasta==0)
				{
					fprintf(fp1,"+%s%s",name1.c_str(),aux1.c_str());//cout<<"+"<<name1<<aux1;
					for (std::vector<char>::iterator it = quality1.begin()+bouderLeft1; it != quality1.end()-bouderRight1; ++it)
    						fprintf(fp1, "%c", *it);//std::cout << *it;
					fprintf(fp1, "\n");//cout<<endl;
					//print quality as integer values 
					//for (std::vector<unsigned int>::iterator it = quality1_int.begin()+bouderLeft1; it != quality1_int.end()-bouderRight1; ++it)
    						//fprintf(fp1, "%d", *it);//std::cout << *it;
					//fprintf(fp1, "\n");//cout<<endl;
				}
				

				name2.erase (0,1);
				//print name
				fprintf(fp2,"%c%s%s",firstC,name2.c_str(),aux2.c_str());//cout<<firstC<<name<<aux1;
				//print bases
				for (std::vector<char>::iterator it = bases2.begin()+bouderLeft2; it != bases2.end()-bouderRight2; ++it)
    					fprintf(fp2, "%c", *it);//std::cout << *it;
				fprintf(fp2, "\n");//cout<<endl;
				//print qualities
				if(args->enableFasta==0)
				{
					fprintf(fp2,"+%s%s",name2.c_str(),aux2.c_str());//cout<<"+"<<name1<<aux1;
					for (std::vector<char>::iterator it = quality2.begin()+bouderLeft2; it != quality2.end()-bouderRight2; ++it)
    						fprintf(fp2, "%c", *it);//std::cout << *it;
					fprintf(fp2, "\n");//cout<<endl;
					//print quality as integer values 
					//for (std::vector<unsigned int>::iterator it = quality2_int.begin()+bouderLeft2; it != quality2_int.end()-bouderRight2; ++it)
    						//fprintf(fp2, "%d", *it);//std::cout << *it;
					//fprintf(fp2, "\n");//cout<<endl;
				}
				
				name1="", name2="", aux1="", aux2=""; 
				bases1.clear();
				bases2.clear();
				quality1.clear();
				quality2.clear();
				quality1_int.clear();
				quality2_int.clear();	
				state=s0;
				break;
		}
	}
	infile1.close();
	infile2.close();
	if ( args->outputFile ) 
	{
		fclose(Outfastq_fp1);
		fclose(Outfastq_fp2);
		fclose(Singleton_fp1);
		fclose(Singleton_fp2);
	}
}

void construct_fastq_trimmed(FILE *fp, char *name, char *bases, uint8_t *quality, int nbases, struct Args * args, struct Names * names)
{

    int j=0;
    uint8_t quality_char;
    int bouderLeft=nbases-1; //por si no encuentra ninguna		
    int bounderRight=nbases-1;

    int sum=0; int k=0; int next=1;
	
    // cut quality values 
    //left
    while (j < nbases && next)
    {
        sum=0;k=0;
	while((j+args->windows_size)<nbases && k<args->windows_size && quality[j+k]>=args->left_cut)
		sum++,k++;

	if(sum==args->windows_size)
		bouderLeft=j, next=0;
	j++;
    }	
    //right
    j=nbases-1; next=1;
    while (j > 0 && next)
    {
	sum=0;k=0;
	while((j-args->windows_size)>0 && k<args->windows_size && quality[j-k]>=args->right_cut)
	   	sum++,k++;
	
	if(sum==args->windows_size)
		bounderRight=j, next=0;
	j--;
    }	

    int readLenght=(bounderRight+1)-bouderLeft;
    if(args->enableFasta==0 && next==0 && readLenght>=args->min_len && readLenght<=args->max_len)
    {
	// print out the name/sequence blocks
    	fprintf(fp, "@%s\n", name);
    
    	for (j = bouderLeft; j <= bounderRight; j++) 
		fprintf(fp,"%c", bases[j]);
    	fprintf(fp, "\n");

    	fprintf(fp, "+%s\n", name);	

    	//print quality as integer values 
    	//for (j = bouderLeft; j <= bounderRight; j++) 
        //	fprintf(fp,"%d ", quality[j] );
    	//fprintf(fp,"\n");
    	

    	//print out quality values (as characters)
    	// formula taken from http://maq.sourceforge.net/fastq.shtml
    	///
    	for (j = bouderLeft; j <= bounderRight; j++) {
        	quality_char = (quality[j] <= 93 ? quality[j] : 93) + 33;
        	fprintf(fp, "%c", (char) quality_char );
    	}	
    	fprintf(fp, "\n");
     }
     else if(next==0 && readLenght>=args->min_len && readLenght<=args->max_len)
     {	
       // print out the name/sequence blocks
    	fprintf(fp, ">%s\n", name);
 	for (j = bouderLeft; j <= bounderRight; j++)
                fprintf(fp,"%c", bases[j]);
        fprintf(fp, "\n");

     }
	
}
void process_sff_to_fastq(struct Args * args, struct Names * names)
{
    sff_common_header h;
    sff_read_header rh;
    sff_read_data rd;
    FILE *sff_fp, *fastq_fp;

    char myArray[names->sffFile.size()+1];
    strcpy(myArray, names->sffFile.c_str());	
    if ( (sff_fp = fopen(&myArray[0], "r")) == NULL )
        cout<<"Could not open file for reading "<<names->sffFile<<endl, exit(1);
    
	
    read_sff_common_header(sff_fp, &h);
    //verify_sff_common_header(PRG_NAME, VERSION, &h);
	
    if ( !args->outputFile ) 
        fastq_fp = stdout;
    
    else 
    {
	if(args->enableFasta==0)	
		names->outputFile+=".fastq";
   	else
		names->outputFile+=".fasta";
	
	myArray[names->outputFile.size()+1];
    	strcpy(myArray, names->outputFile.c_str());
        if ( (fastq_fp = fopen(myArray, "w")) == NULL ) 
        	cout<<"Could not open file for writing "<<names->outputFile<<endl, exit(1);
    }

    int left_clip = 0, right_clip = 0, nbases = 0, trim_flag= 1;
    char *name;
    char *bases;
    uint8_t *quality;
    register int i;
    int numreads = (int) h.nreads;
    for (i = 0; i < numreads; i++) {
        read_sff_read_header(sff_fp, &rh);
        read_sff_read_data(sff_fp, &rd, h.flow_len, rh.nbases);

        // get clipping points 
        get_clip_values(rh, trim_flag, &left_clip, &right_clip);
        nbases = right_clip - left_clip;

        // create bases string 
        bases = get_read_bases(rd, left_clip, right_clip);

        // create quality array 
        quality = get_read_quality_values(rd, left_clip, right_clip);

        // create read name string 
        int name_length = (int) rh.name_len + 1; // account for NULL termination
        name = (char *) malloc( name_length * sizeof(char) );
        if (!name) {
            fprintf(stderr, "Out of memory! For read name string!\n");
            exit(1);
        }
        memset(name, '\0', (size_t) name_length);
        strncpy(name, rh.name, (size_t) rh.name_len);
        
        if (strlen(bases) > 0) 
		construct_fastq_trimmed(fastq_fp, name, bases, quality, nbases, args, names);	
        

        free(name);
        free(bases);
        free(quality);
        free_sff_read_header(&rh);
        free_sff_read_data(&rd);
    }

    free_sff_common_header(&h);
    fclose(fastq_fp);
    fclose(sff_fp);
}

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

        }

	if (error)
		help_message(), exit(0);
}

void help_message() 
{
    cout<< "Usage: "<<PRG_NAME<<" [options]"<<endl;//"[454 sff_file or illumina fastq_files]\n";
    cout<< "-h"<<"\t\t\tThis help message\n";
    cout<< "-v"<<"\t\t\tProgram and version information\n";
    cout<< "-f"<<"\t\t\tEnable fasta output (default fastq)\n";
    cout<< "-l<int>"<<"\t\t\tSet lef-cut value for quality scores (int default 0)\n";
    cout<< "-r<int>"<<"\t\t\tSet right-cut value for quality scores (int default 0)\n";
    cout<< "-w<int>"<<"\t\t\tSet windows size to check on the quality scores (int default 1)\n";
    cout<< "-m<int>"<<"\t\t\tFilter sequence shorter than min_len (int default 1)\n";
    cout<< "-mx<int>"<<"\t\tFilter sequence larger than max_len (int default 5000)\n";
    cout<< "-sff<454 files>"<<"\t\t454 input (file.sff)\n";
    cout<< "-i<illumina files>"<<"\tIllumina inputs(file1.fastq file2.fastq)\n";
    cout<< "-o<fastq_file>"<<"\t\tDesired fastq output file."<<" If not specified to stdout\n";
    cout<< "Quality input and output is calculated according the formula taken from http://maq.sourceforge.net/fastq.shtml"<<endl;	
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

