//
//  gRNA.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#ifndef __Casper_Seq_Finder__gRNA__
#define __Casper_Seq_Finder__gRNA__

#include <iostream>
#include <vector>
#include <set>
using namespace std;


class gRNA {
public:
    gRNA();
    ~gRNA();
    
public:
    unsigned long insertSequence(long, int, int, bool, bool, std::string, int);  //returns the seed sequence to the map after processing
    long getLocation() {return PAMlocation;} // the negative refers to the strand direction +/-: sense/anti-sense
    std::string getHypLoc();
    std::string getHypTail();
    std::string getHypPam();
    std::pair<unsigned long, std::string> getVectorPair(unsigned long,bool);
    int chrNumber() {return Chromosome;}
    
    std::string baseConvert(unsigned long long, int);
    
    std::string getScore() {return baseConvert(OnScore,64);}

    
private:
    unsigned int tailSeq; // capable of storing a 8 nucleotide sequence (probably will only use half)
    unsigned int pamSeq; // capable of storing an 8 nucleotide PAM sequence
    long PAMlocation;
    int Chromosome;
    int Pamsize;
    bool Anti;
    int SeedLength;
    int OnScore;
    
private:
    std::string decompressSeq(unsigned long, int);
    unsigned long compressSeq(std::string);
    int convertCharBase4(char);
    char convertBase4toChar(int);
    
};


#endif /* defined(__Casper_Seq_Finder__gRNA__) */
