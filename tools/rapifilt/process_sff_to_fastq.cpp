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
#include "sff.h"
#include "process_sff_to_fastq.h"

void construct_fastq_trimmed(FILE *fp, char *name, char *bases, unsigned int *quality, int nbases, struct Args * args, struct Names * names, int *bounders)
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
            //fprintf(fp,"%d ", quality[j] );
        //fprintf(fp,"\n");

        //print out quality values (as characters)
        // formula taken from http://maq.sourceforge.net/fastq.shtml
        for (j = bouderLeft; j <= bounderRight; j++) 
        {
            quality_char = (quality[j] <= 93 ? quality[j] : 93) + 33;
            fprintf(fp, "%c", (char) quality_char );
        }
        fprintf(fp, "\n");
     }
     else if(next==0 && readLenght>=args->min_len && readLenght<=args->max_len)
     {
        //print out the name/sequence blocks
        fprintf(fp, ">%s\n", name);
        for (j = bouderLeft; j <= bounderRight; j++)
            fprintf(fp,"%c", bases[j]);
        fprintf(fp, "\n");

     }
     
    bounders[0]=bouderLeft;
    bounders[1]=bounderRight;

}

void process_sff_to_fastq(struct Args * args, struct Names * names)
{
    sff_common_header h;
    sff_read_header rh;
    sff_read_data rd;
    FILE *sff_fp, *fastq_fp, *stat_fp;

    char myArray[names->sffFile.size()+1], myArray2[names->sffFile.size()+1];
    strcpy(myArray, names->sffFile.c_str());	
    if ( (sff_fp = fopen(&myArray[0], "r")) == NULL )
        cout<<"Could not open file for reading "<<names->sffFile<<endl, exit(1);
    

    read_sff_common_header(sff_fp, &h);
    //verify_sff_common_header(PRG_NAME, VERSION, &h);

    if ( !args->outputFile ) 
        fastq_fp = stdout, stat_fp=stdout;
    
    else 
    {
        myArray2[(names->outputFile+"_stats.txt").size()+1];
        strcpy(myArray2, (names->outputFile+"_stats.txt").c_str());
        
        if(args->enableFasta==0)	
            names->outputFile+=".fastq";
        else
            names->outputFile+=".fasta";

        myArray[names->outputFile.size()+1];
        strcpy(myArray, names->outputFile.c_str());
        
        
        if ( (fastq_fp = fopen(myArray, "w")) == NULL  || (stat_fp = fopen(myArray2, "w")) == NULL) 
            cout<<"Could not open file for writing "<<names->outputFile<<endl, exit(1);
    }

    int left_clip = 0, right_clip = 0, nbases = 0, trim_flag= 1;
    char *name;
    char *bases;
    unsigned int *quality;
    register int i;
    int numreads = (int) h.nreads;
    
    //statistics
    int bounders [2]={};
    int binSize=args->binSize;
    int Gbase=0, Cbase=0;
    int Gbasef=0, Cbasef=0;
    int totalBases=0, totalBasesf=0;
    int qual_position[maxReadLen]={0};  //maxReadLen max read lenght
    int qual_positionf[maxReadLen]={0};  //maxReadLen max read lenght
    int maxQ[maxReadLen]={0}, minQ[maxReadLen];
    int maxQF[maxReadLen]={0}, minQF[maxReadLen];   
    int total_reads[maxReadLen]={0}, total_readsf[maxReadLen]={0};   
    int max_read_len=0, max_read_lenf=0;

    
    for (i = 0; i < numreads; i++) 
    {
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
        {
            //compute orginal statistics
            statistic_mean (quality, binSize,nbases, &max_read_len, &maxQ[0], &minQ[0], &total_reads[0], &qual_position[0],0,0);
            //trimer
            construct_fastq_trimmed(fastq_fp, name, bases, quality, nbases, args, names,&bounders[0]);
            //final statistic_mean
            statistic_mean (quality, binSize,nbases-bounders[1], &max_read_lenf, &maxQF[0], &minQF[0], &total_readsf[0], &qual_positionf[0],bounders[0],bounders[1]);

        }

        free(name);
        free(bases);
        free(quality);
        free_sff_read_header(&rh);
        free_sff_read_data(&rd);
    }
    
    //statistics process  for original reads
    fprintf(stat_fp,">Original Reads: %d%s%d%s", total_reads[0],"\tLargest read: ",max_read_len,"\n");
    int index=max_read_len/binSize;
    for(int i=0;i<index;i++)
        fprintf(stat_fp, "%d%s", qual_position[i]/total_reads[i]," ");//std::cout << *it;
    fprintf(stat_fp, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp, "%d%s", maxQ[i]," ");//std::cout << *it;
    fprintf(stat_fp, "\n");//cout<<endl;

    for(int i=0;i<index;i++)
        fprintf(stat_fp, "%d%s", minQ[i]," ");//std::cout << *it;                     
    fprintf(stat_fp, "\n");//cout<<endl;
    
    //statistics process  for final reads
    fprintf(stat_fp,">Final Reads: %d%s%d%s", total_readsf[0],"\tLargest read: ",max_read_lenf,"\n");
    int indexF=max_read_lenf/binSize;
    for(int i=0;i<indexF;i++)
        fprintf(stat_fp, "%d%s", qual_positionf[i]/total_readsf[i]," ");//std::cout << *it;
    fprintf(stat_fp, "\n");//cout<<endl;

    for(int i=0;i<indexF;i++)
        fprintf(stat_fp, "%d%s", maxQF[i]," ");//std::cout << *it;
    fprintf(stat_fp, "\n");//cout<<endl;

    for(int i=0;i<indexF;i++)
        fprintf(stat_fp, "%d%s", minQF[i]," ");//std::cout << *it;                     
    fprintf(stat_fp, "\n");//cout<<endl;
    
    free_sff_common_header(&h);
    fclose(fastq_fp);
    fclose(sff_fp);
}






