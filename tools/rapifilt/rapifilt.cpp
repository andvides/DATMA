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

#include "rapifilt.h"


/* M A I N *******************************************************************/
int main(int argc, char *argv[]) 
{
    Args   args = {false, false, false, false, false, false, 0, 0, 1, 1,5000,0,0,1};
    Names  names;

    process_options(argc, argv, &args, &names);
    

    if(args.illuminaFile)
        process_illumina_to_fastq(&args, &names);
    else if(args.sffFile)	
    	process_sff_to_fastq(&args, &names);
    else if(args.gz)	
    	process_gzFile(&args, &names);
    else if(args.fastq)	
    	process_fastqFile(&args, &names);
    else
        help_message(),	exit(0);
    return 0;
}

