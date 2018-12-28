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

#include "process_gzFile.h"
#include "generic.h"

using namespace std;

void process_gzFile(struct Args * args, struct Names * names)
{
    std::ifstream file(names->gz_illuminaFile1.c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    //Convert streambuf to istream
    std::istream instream(&inbuf);
    //Iterate lines
    std::string line;
    while(std::getline(instream, line)) {
        std::cout << line << std::endl;
    }
    //Cleanup
    file.close();
}