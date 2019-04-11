#include<stdio.h> 
#include <stdlib.h>
#define MAXCHAR 1000
//#define joinWith "??????????????"
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
int matchNames(char * str1,  char *str2)
{
    //printf("matching name of files\n");
    printf("%s",str1);
    while(*str1++ != ' ' && *str2++ != ' ') //stop in the first blank space
    {
        if(*str1 != *str2)
        {
             printf("f1 and r1 must have same name\n");
             return 10; //error            
        }
    }
    return 1; //ok
}
int extendedReads(char * str1,  char *str2, char *joinString)
{
    //printf("Extending\n");
    int bases=strlen(str1);
    int b=0;
    for(b=0;b<bases-1;b++) //to not take the enter
            printf("%c",str1[b]);
    printf("%s",joinString);
    
    bases=strlen(str2);
    char reverse;
    for(b=bases-2;b>=0;b--)//NOt take the enter
    {
        reverse=nt_table_reverse[str2[b]];
        printf("%c",reverse);
    }
    printf("\n");
    
    return 2; //ok
}
int extendedQualities(char * str1,  char *str2, char *joinString)
{
    //printf("Extending\n");
    int bases=strlen(str1);
    int b=0;
    for(b=0;b<bases-1;b++) //to not take the enter
        printf("%c",str1[b]);
    
    bases=strlen(joinString);
    for(b=0;b<bases;b++) //to take the enter
        printf("%c",'!');//minimal quality

    bases=strlen(str2);
    for(b=bases-2;b>=0;b--)//NOt take the enter
        printf("%c",str2[b]);
    printf("\n");
    return 0; //ok
}
void loadFiles(FILE *f1, FILE *r1, char *joinString)
{
    char str1[MAXCHAR];
    char str2[MAXCHAR];
    int state=0;
    //printf("Reading files\n");
    while (fgets(str1, MAXCHAR, f1) != NULL)
    {   
        
        if(fgets(str2, MAXCHAR, r1) == NULL)
            printf("f1 and r1 do not have the same number of lines\n"), exit(2);
        
       if (state==0) //read first character
            state=matchNames(&str1[0],&str2[0]);
        else if (state==1) //read bases
            state=extendedReads(&str1[0],&str2[0],joinString);//+
        else if (state==2)
             printf("%s",str1),state=3;
        else if (state==3)// read qualities
            state=extendedQualities(&str1[0],&str2[0],joinString);
        else
           printf("Error reading files\n"),exit(3);
    }
        
    fclose(f1); 
    fclose(r1); 
}

int main(int argc,char* argv[]) 
{ 
    int counter; 
    if(argc<3 || argc > 4) 
        printf("Use: mergetNotCombined F1.fastq R1.fastq joinString\n"), exit(1);
    
    
    FILE *f1 = fopen(argv[1],"r"); 
    FILE *r1 = fopen(argv[2],"r"); 
    
    if(f1 == NULL || r1==NULL)
        printf("File doesn´t exits\n"),  exit(1);
    
    loadFiles(f1,r1,argv[3]);
    
    return 0; 
}
