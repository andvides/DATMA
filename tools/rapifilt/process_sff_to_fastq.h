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


void construct_fastq_trimmed(FILE *fp, char *name, char *bases, uint8_t *quality, int nbases, struct Args * args, struct Names * names, int *bounders);
void process_sff_to_fastq(struct Args * args, struct Names * names);

