//
//  CrisprGroup.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#ifndef __Casper_Seq_Finder__CrisprGroup__
#define __Casper_Seq_Finder__CrisprGroup__

#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include "hashset.h"
#include "hashmap.h"
#include "gRNA.h"
#include "set.h" //belongs to the StanfordCPPLibrary Copyright (c) Stanford University 2008.
#include <map>
#include <unordered_map>
#include "foreach.h"
#include "pameval.h"
#include "Scoring.h"

class CrisprGroup {
public:
    CrisprGroup(int, std::string, std::string);
    ~CrisprGroup();
public:
    void findPAMs(std::string &s, bool, int, std::string pam, bool, bool, std::string scr);
    void printSequences();
    void processTargets();  //CAN ONLY BE CALLED ONCE AND AT THE END OF THE SEARCH PROCESS!!!
private:
    std::string filename;
private:
    void addToMap(unsigned long, gRNA*);
    int charToInt(char);
    void tempFileNAGs(int);
    static pamEval evaluate(std::string, bool);
    std::string compressSequence(std::string);
    gRNA* sCur;
    int numChromosomes;
    int curchrom;
    
private:
    std::unordered_map<unsigned long, std::vector<gRNA*>> Seed_Map; //Stores all the potential target sites
    std::vector<std::vector<std::pair<long, std::string>>> total_seqs; //sorted unique sequences
    std::vector<std::pair<unsigned long, std::vector<gRNA*>>> repeat_seqs; //unsorted repeated sequences
    
/* Stuff for iteration */
public:
    int chrCount() {return numChromosomes;}
    unsigned long Size(int k) {return total_seqs.at(k).size();}
    unsigned long totSize();
    unsigned long repSize() {return repeat_seqs.size();}
    std::string nextUnique(int, long);
    std::pair<unsigned long, std::vector<gRNA*>> nextRepeatSet(int i) {return repeat_seqs[i];}
};


#endif /* defined(__Casper_Seq_Finder__CrisprGroup__) */
