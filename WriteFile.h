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
#include <ctype.h>
#include "CrisprGroup.h"
#include "gRNA.h"
#include "pameval.h"
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "sqlite3.h"

using namespace std;

class WriteFile {
public:
    WriteFile();
    ~WriteFile();
public:
    void inputStats(std::vector<int>,std::string);
    void setFileName(string, string);
    void retrieveData(CrisprGroup*,std::vector<std::string>, bool repeats);
private:
    void inputRepeatData(gRNA* g);
    void entry();
    int charToInt(char c);
    char convertBase4toChar(int);
    std::string decompressSeq(unsigned long long, short);
private:
    pamEval PAMstat;
    string filename;
    string chr_stats_str;
    string mystats;
    ofstream outputfile;
    ofstream repeatfile;
    int chromosome;
    std::string position;
    char strand;
    string sequence;
    string score;
    vector<int> chromosomeseqcount;
};


#endif /* defined(__Casper_Seq_Finder__WriteFile__) */
