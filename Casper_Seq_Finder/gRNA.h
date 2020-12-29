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
    unsigned long long insertSequence(long, int, int, bool, std::string, int, short, short, short);  //returns the seed sequence to the map after processing
    long getLocation() {return PAMlocation;} // the negative refers to the strand direction +/-: sense/anti-sense
    unsigned int getFiveSeq() {return fiveSeq;}
    unsigned int getThreeSeq() {return threeSeq;}
    unsigned short getPam() {return pamSeq;}
    int chrNumber() {return Chromosome;}
    short getScore() {return OnScore;}

    
private:
    unsigned long long compressSeq(std::string);
    int convertCharBase4(char);
    unsigned int fiveSeq; // capable of storing a 16 nucleotide sequence
    unsigned int threeSeq; // capable of storing a 16 nucleotide sequence
    unsigned short pamSeq; // capable of storing an 8 nucleotide PAM sequence
    long PAMlocation;
    unsigned short Chromosome;
    short OnScore;
    
};


#endif /* defined(__Casper_Seq_Finder__gRNA__) */
