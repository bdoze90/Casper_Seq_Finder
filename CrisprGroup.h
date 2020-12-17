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
#include "gRNA.h"
#include <map>
#include <unordered_map>
#include "pameval.h"
#include "Scoring.h"

class CrisprGroup {
public:
    CrisprGroup(int, std::string, std::string,int,int);
    ~CrisprGroup();
public:
    void findPAMs(std::string &s, bool, int, std::string pam, bool, bool, std::string scr);
    void printSequences();
    void processTargets();  //CAN ONLY BE CALLED ONCE AND AT THE END OF THE SEARCH PROCESS!!!
	int pam_length;
	int len_seed;
	int len_seq;
private:
    std::string filename;
    
private:
    void addToMap(unsigned int, gRNA*);
    int charToInt(char);
    static pamEval evaluate(std::string, bool);
    std::string compressSequence(std::string);
	char convertBase4toChar(int i);
    gRNA* sCur;
    int numChromosomes;
    int curchrom;
    
private:
    std::unordered_map<unsigned int, std::vector<gRNA*>> Seed_Map; //Stores all the potential target sites
    std::vector<std::vector<std::pair<long, std::string>>> total_seqs; //sorted unique sequences
    std::vector<std::pair<unsigned int, std::vector<gRNA*>>> repeat_seqs; //unsorted repeated sequences
    
/* Stuff for iteration */
public:
    int chrCount() {return numChromosomes;}
    unsigned long Size(int k) {return total_seqs.at(k).size();}
    unsigned long totSize();
    unsigned long repSize() {return repeat_seqs.size();}
    std::string nextUnique(int, long);
    std::pair<unsigned int, std::vector<gRNA*>> nextRepeatSet(int i) {return repeat_seqs[i];}
	string decompress64(string seq, int length);
	std::string decompressSeq(unsigned long cseq, short exp_len);
	unsigned int decompress_ontarget(string base64seq);
	string int2nt(int num);
};


#endif /* defined(__Casper_Seq_Finder__CrisprGroup__) */
