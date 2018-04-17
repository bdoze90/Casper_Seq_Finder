//
//  WriteFile.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#ifndef __Casper_Seq_Finder__WriteFile__
#define __Casper_Seq_Finder__WriteFile__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "CrisprGroup.h"
#include "gRNA.h"

using namespace std;

class WriteFile {
public:
    WriteFile();
    ~WriteFile();
public:
    void setFileName(string fn);
    void retrieveData(CrisprGroup*);
    void printInfo(CrisprGroup*);
private:
    void inputData(gRNA* g);
    void entry();
    int charToInt(char c);
    std::string compressSequence(std::string s);
private:
    string filename;
    ofstream outputfile;
    int chromosome;
    std::string position;
    char strand;
    string sequence;
    string score;
    vector<int> chromosomeseqcount;
};


#endif /* defined(__Casper_Seq_Finder__WriteFile__) */
