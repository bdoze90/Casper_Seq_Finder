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
#include "RepStor.h"

class CrisprGroup {
public:
    CrisprGroup(unsigned long, pamEval, std::string,std::string, bool);
    ~CrisprGroup();
public:
    void findPAMs(std::string &s, bool, int, std::string scr);
    void findPAMs_notRepeats(std::string &s, bool, int, std::string scr);
    void printSequences();
    void initiateTotalSeqs();
    void processTargets(std::vector<int>);  //CAN ONLY BE CALLED ONCE AND AT THE END OF THE SEARCH PROCESS!!!
private:
    std::string filename;
private:
    void addToMap(unsigned long, gRNA*,bool);
    int charToInt(char);
    static pamEval evaluate(std::string, bool);
    std::string compressSequence(std::string);
    gRNA* sCur;
    int numChromosomes;
    int curchrom;
    pamEval PAMstat;
    struct compgrna {
        unsigned long seed;
        unsigned int cfive;
        unsigned int cthree;
        short cpam;
        short score;
        bool strand = true;
    };
private:  // Functions for manipulating the output of the gRNA objects
    static bool pairCompare(const std::pair<long, compgrna>, const std::pair<long, compgrna>);
    std::string decompressSeq(unsigned long, short);
    char convertBase4toChar(int);
    std::string baseConvert(unsigned long long, int);
    
private:
    std::unordered_map<unsigned long, std::vector<gRNA*>> Seed_Map; //Stores all the potential target sites
    std::vector<std::vector<std::pair<long, compgrna>>> total_seqs; //sorted unique sequences
    std::vector<std::pair<unsigned int, std::vector<gRNA*>>> repeat_seqs; //unsorted repeated sequences
    
/* Stuff for iteration */
public:
    int chrCount() {return numChromosomes;}
    pamEval getPamEval() {return PAMstat;}
    unsigned long Size(int k) {return total_seqs.at(k).size();}
    unsigned long totSize();
    unsigned long repSize() {return repeat_seqs.size();}
    std::pair<long,std::string> nextUnique(int, long);
    std::pair<unsigned int, std::vector<gRNA*>> nextRepeatSet(int i) {return repeat_seqs[i];}
};


#endif /* defined(__Casper_Seq_Finder__CrisprGroup__) */
