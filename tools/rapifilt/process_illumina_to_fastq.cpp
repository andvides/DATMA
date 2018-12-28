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
#include "process_illumina_to_fastq.h"

void process_illumina_to_fastq(struct Args * args, struct Names * names)
{
    
    //output file
    FILE *Outfastq_fp1, *Outfastq_fp2, *Singleton_fp1, *Singleton_fp2, *bad_fp1, *bad_fp2, *fp1, *fp2, *stat_fp1, *stat_fp2;
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
        string stat1=names->outputFile+"_1stat.txt";
        string stat2=names->outputFile+"_2stat.txt";

        char myArray1[name5.size()+1];
        char myArray2[name3.size()+1];		
        char myArray3[single1.size()+1];
        char myArray4[single2.size()+1];
        char myArray5[bad1.size()+1];
        char myArray6[bad2.size()+1];
        char myArray7[stat1.size()+1];
        char myArray8[stat2.size()+1];

        strcpy(myArray1, name5.c_str());
        strcpy(myArray2, name3.c_str());
        strcpy(myArray3, single1.c_str());
        strcpy(myArray4, single2.c_str());
        strcpy(myArray5, bad1.c_str());
        strcpy(myArray6, bad2.c_str());
        strcpy(myArray7, stat1.c_str());
        strcpy(myArray8, stat2.c_str());
        
        if  ((Outfastq_fp1 	= fopen(myArray1, "w")) == NULL || (Outfastq_fp2 	= fopen(myArray2, "w"))	== NULL ||
            (Singleton_fp1 	= fopen(myArray3, "w")) == NULL || (Singleton_fp2 	= fopen(myArray4, "w"))	== NULL ||
            (bad_fp1 		= fopen(myArray5, "w")) == NULL || (bad_fp2 		= fopen(myArray6, "w"))	== NULL ||
            (stat_fp1 		= fopen(myArray7, "w")) == NULL || (stat_fp2	    = fopen(myArray8, "w"))	== NULL) 
                cout<<"Could not open files for writing "<<names->outputFile<<endl, exit(1);
    }

    char * file1 = strdup(names->illuminaFile1.c_str());
    char * file2 = strdup(names->illuminaFile2.c_str());
    //cout<<"Input Files: "<<names->illuminaFile1<<", "<<names->illuminaFile2<<endl;

    //read File
    ifstream infile1, infile2;
    infile1.open(file1);
    infile2.open(file2);

    typedef enum {s0, s1,s2, s3, s4,s5,s6,s7,s8,s9,s10, s11,s12,s13} STATES; //to read multifasta files
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
    
    //statistics vector
    int binSize=args->binSize;
    int Gbase_F=0, Gbase_R=0;
    int Cbase_F=0, Cbase_R=0;
    int Gbase_Ff=0, Gbase_Rf=0;
    int Cbase_Ff=0, Cbase_Rf=0;

    int total_readsF[maxReadLen]={0}, total_readsR[maxReadLen]={0};   
    int total_readsFf[maxReadLen]={0}, total_readsRf[maxReadLen]={0};   

   
    int qual_position_F[maxReadLen]={0};  //maxReadLen max read lenght
    int qual_position_R [maxReadLen]={0}; //maxReadLen max read lenght
    int qual_position_Ff[maxReadLen]={0};  //maxReadLen max read lenght
    int qual_position_Rf [maxReadLen]={0}; //maxReadLen max read lenght
    

    int maxQF[maxReadLen]={0}, minQF[maxReadLen];
    int maxQR[maxReadLen]={0}, minQR[maxReadLen];   
    int maxQFf[maxReadLen]={0}, minQFf[maxReadLen];
    int maxQRf[maxReadLen]={0}, minQRf[maxReadLen];   

    int i=0, max_read_lenF=0, max_read_lenR=0, max_read_lenFf=0,max_read_lenRf=0;
    
    for (int j=0; j<maxReadLen;j++)
        minQF[j]=1000, minQR[j]=1000, minQFf[j]=1000, minQRf[j]=1000;
     
    while (next)
    {
        char word1='\0', word2='\0'; 
        //cout<<state<<endl;
        switch ( state ) 
        {
            case s0:
                next=false;//total_reads++;
                while(word1!='\n') //forward read name
                {
                    if(infile1.get(word1))
                    {
                        next=true;
                        if(int(word1)!=32 && save) 
                            name1+=word1;
                        else
                            aux1+=word1, save=false;
                    }
                    else
                        break;
                }
                state=s1;
                save=true;
                break;
            case s1:
                next=false;
                while(word2!='\n') //Reverse read name
                {	
                    if(infile2.get(word2))
                    {
                        next=true;
                        if(int(word2)!=32 && save) 
                            name2+=word2;
                        else
                            aux2+=word2, save=false;
                    }
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
                while(infile1.get(word1) && word1!='\n') //forward bases
                    bases1.push_back(word1);
                state=s4;
                break;
            case s4:
                while(infile2.get(word2) && word2!='\n') //reverse bases
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
                while(infile1.get(word1) && word1!='\n') //Forward quality
                    quality1.push_back(word1),quality1_int.push_back((unsigned int)word1-33);
                state=s7;
                break;
            case s7:
                while(infile2.get(word2) && word2!='\n') //Reverse quality
                    quality2.push_back(word2),quality2_int.push_back((unsigned int)word2-33) ;;
                state=s8;
                break;
            case s8:
                //cut by quality forward
                
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
                //cut by quality reverse
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
                        fp1=Outfastq_fp1, fp2=Outfastq_fp2;//,total_readsf++;
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
                    //fprintf(fp1,"+%s%s",name1.c_str(),aux1.c_str());//cout<<"+"<<name1<<aux1;
                    fprintf(fp1,"+\n");//cout<<"+"<<name1<<aux1;
                    for (std::vector<char>::iterator it = quality1.begin()+bouderLeft1; it != quality1.end()-bouderRight1; ++it)
                        fprintf(fp1, "%c", *it);//std::cout << *it;
                    fprintf(fp1, "\n");//cout<<endl;
                    //print quality as integer values 
                    //for (std::vector<unsigned int>::iterator it = quality1_int.begin()+bouderLeft1; it != quality1_int.end()-bouderRight1; ++it)
                        //fprintf(fp1, "%d%s", *it," ");//std::cout << *it;
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
                    //fprintf(fp2,"+%s%s",name2.c_str(),aux2.c_str());//cout<<"+"<<name1<<aux1;
                    fprintf(fp2,"+\n");
                    for (std::vector<char>::iterator it = quality2.begin()+bouderLeft2; it != quality2.end()-bouderRight2; ++it)
                        fprintf(fp2, "%c", *it);//std::cout << *it;
                    fprintf(fp2, "\n");//cout<<endl;
                    //print quality as integer values 
                    //for (std::vector<unsigned int>::iterator it = quality2_int.begin()+bouderLeft2; it != quality2_int.end()-bouderRight2; ++it)
                        //fprintf(fp2, "%d%s", *it," ");//std::cout << *it;
                    //fprintf(fp2, "\n");//cout<<endl;
                }
                
                state=s12;
                break;
            case s12: //statistics process  for original read
                
                 //original
                GC_content(bases1.data(),size1,&Gbase_F,&Cbase_F,0,0); //forward
                GC_content(bases2.data(),size1,&Gbase_R,&Cbase_R,0,0); // reverse
                
                statistic_mean (quality1_int.data(), binSize, size1, &max_read_lenF, &maxQF[0], &minQF[0], &total_readsF[0], &qual_position_F[0], 0, 0);
                statistic_mean (quality2_int.data(), binSize, size1, &max_read_lenR, &maxQR[0], &minQR[0], &total_readsR[0], &qual_position_R[0], 0, 0);

                state=s13;
                break;
            case s13: //statistics process  for trimmer read
                
                //Final
                GC_content(bases1.data(),size1,&Gbase_Ff,&Cbase_Ff,bouderLeft1,bouderRight1); //forward
                GC_content(bases2.data(),size1,&Gbase_Rf,&Cbase_Rf,bouderLeft2,bouderRight2); // reverse
                
                statistic_mean (quality1_int.data(), binSize, finalSize1, &max_read_lenFf, &maxQFf[0], &minQFf[0], &total_readsFf[0], &qual_position_Ff[0], bouderLeft1, bouderRight1);
                statistic_mean (quality2_int.data(), binSize, finalSize2, &max_read_lenRf, &maxQRf[0], &minQRf[0], &total_readsRf[0], &qual_position_Rf[0], bouderLeft2, bouderRight2);

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
    
    
    //statistics process  for original reads
    fprintf(stat_fp1,">Original Reads: %d%s%d%s", total_readsF[0],"\tLargest read: ",max_read_lenF,"\n");
    int index=max_read_lenF/binSize;
    for(int i=0;i<index;i++)
        fprintf(stat_fp1, "%d%s", qual_position_F[i]/total_readsF[i]," ");//std::cout << *it;
    fprintf(stat_fp1, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp1, "%d%s", maxQF[i]," ");//std::cout << *it;
    fprintf(stat_fp1, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp1, "%d%s", minQF[i]," ");//std::cout << *it;                     
    fprintf(stat_fp1, "\n");//cout<<endl;
    
    fprintf(stat_fp2,">Original Reads: %d%s%d%s", total_readsR[0],"\tLargest read: ",max_read_lenR,"\n");
    index=max_read_lenR/binSize;
    for(int i=0;i<index;i++)
        fprintf(stat_fp2, "%d%s", qual_position_R[i]/total_readsR[i]," ");//std::cout << *it;
    fprintf(stat_fp2, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp2, "%d%s", maxQR[i]," ");//std::cout << *it;
    fprintf(stat_fp2, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp2, "%d%s", minQR[i]," ");//std::cout << *it;                     
    fprintf(stat_fp2, "\n");//cout<<endl;
    
    //statistics process  for final reads
    fprintf(stat_fp1,">Final Reads: %d%s%d%s", total_readsFf[0],"\tLargest read: ",max_read_lenFf,"\n");
    index=max_read_lenFf/binSize;
    for(int i=0;i<index;i++)
        fprintf(stat_fp1, "%d%s", qual_position_Ff[i]/total_readsFf[i]," ");//std::cout << *it;
    fprintf(stat_fp1, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp1, "%d%s", maxQFf[i]," ");//std::cout << *it;
    fprintf(stat_fp1, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp1, "%d%s", minQFf[i]," ");//std::cout << *it;                     
    fprintf(stat_fp1, "\n");//cout<<endl;
    
    fprintf(stat_fp2,">Original Reads: %d%s%d%s", total_readsRf[0],"\tLargest read: ",max_read_lenRf,"\n");
    index=max_read_lenRf/binSize;
    for(int i=0;i<index;i++)
        fprintf(stat_fp2, "%d%s", qual_position_Rf[i]/total_readsRf[i]," ");//std::cout << *it;
    fprintf(stat_fp2, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp2, "%d%s", maxQRf[i]," ");//std::cout << *it;
    fprintf(stat_fp2, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp2, "%d%s", minQRf[i]," ");//std::cout << *it;                     
    fprintf(stat_fp2, "\n");//cout<<endl;
    
    
    
    infile1.close();
    infile2.close();
    if ( args->outputFile ) 
    {
        fclose(Outfastq_fp1);
        fclose(Outfastq_fp2);
        fclose(Singleton_fp1);
        fclose(Singleton_fp2);
        fclose(stat_fp1);
        fclose(stat_fp2);
    }
}

