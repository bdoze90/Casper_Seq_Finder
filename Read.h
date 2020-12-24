//
//  Read.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#ifndef __Casper_Seq_Finder__Read__
#define __Casper_Seq_Finder__Read__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

class Read {
public:
    Read() { }
    void setFileName(std::string name) { filename = name; }
    void openFile();
    std::string getLine();
    std::string FirstLine();
    bool newLine();
    void closeFile();
    //For reading the finder code file
public:
    std::string getPAM();
    std::vector<std::string> getOPAMs();
    std::string getOrgCode();
    std::vector<std::string> getFileLocations();
    bool getAnti();
    int getOpamNum() { return opamnum; }
private:
    int opamnum;
    std::string filename;
    std::ifstream* stream;
};


#endif /* defined(__Casper_Seq_Finder__Read__) */

