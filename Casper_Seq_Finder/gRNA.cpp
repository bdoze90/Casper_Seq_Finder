//
//  gRNA.cpp
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include "gRNA.h"
#include <vector>
#include <cmath>


/* Constructor: gRNA
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Initializes a gRNA object setting the appropriate parts of the object to their default value
 */

gRNA::gRNA() {
    PAMlocation = 0;
    Chromosome = 0;
}

/* Destructor: gRNA
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Orderly deletes information stored in the gRNA object by careful iteration through pointers and
 * associated linked list in the form of the RelativeList (calls another destructor).
 */

gRNA::~gRNA() {
    PAMlocation = 0;
    Chromosome = 0;
    Pamsize = 0;
    tailSeq = NULL;
    //does there need to be any extra here i.e. deleting pointers?
}


/* Function: insertSequence
 * ---------------------------------------------------------------------------------------------------------
 * Usage: After an instance of gRNA is created, this function assigns the crucial parameters to the internal
 * assignments, being the complete sequence (g_sequence) and the location of the sequence stored as the
 * location of the protospacer adjacent motif (PAMlocation). This function also returns the seed sequence so
 * that it can be stored in the Seed_Map of CrisprGroup for comparison to discover repeats.
 */

unsigned long gRNA::insertSequence(long index, int chr, int pamsize, bool anti, bool strand, std::string seq, int score) {
    PAMlocation = index;
    if (!strand) {
        PAMlocation *= -1;
    }
    OnScore = score;
    Chromosome = chr+1; // this added 1 takes care of the indexing error and properly assigns chromosome number.
    Pamsize = pamsize;
    Anti = anti;  // !this needs to be judged when inputting tail and seed seqs unless taken care of somewhere else
    /* The following code process the total sequence into multiple elements including the tailSeq, seedSeq, etc. */
    std::string seed;
    std::string tail;
    seed = seq.substr(seq.size()-16,16);  // 16 is the maximum number of nucleotides able to be stored in unsigned long seed.
    tail = seq.substr(0,seq.size()-16);
    tailSeq = compressSeq(tail); //should be fine because limit is 8 nucleotides
    return compressSeq(seed);
}

/* Function: pamSequence  !!!WRONG!!!
 * ---------------------------------------------------------------------------------------------------------
 * Usage: returns the part of the total sequence sequence that is the PAM.
 */

std::string gRNA::pamSequence(std::string fullseq) {
    std::string pam;
    if (Anti) {
        pam = fullseq.substr(0,Pamsize);
    } else {
        pam = fullseq.substr(fullseq.size()-Pamsize,20);
    }
    return pam;
}

/* Function: decompressSeq
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Takes the stored sequence (which is compressed to base-10) and outputs it as its string of A,T,C,G.
 * There is a problem with an Adenine in the first part of the seedSeq or tailSeq where it is not taken into
 * account because of the fact that it is represented by a 0.  To remedy this, when the decompression occurs
 * we need to make sure that if the decompression reveals a sequence shorter than the expected length, to add
 * an A to the end.
 */

std::string gRNA::decompressSeq() {
    std::string uncompressed;
    //do the reverse binary transition from base-10 to base-4
    return uncompressed;
}

/* Function: compressSeq
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Called from the insertSequence function to take in the passed sequence and compress it into a
 * base-10 integer that can be stored with data requirements.
 * Special note: THIS PROCESS FLIPS THE STRING DIRECTION SO A SEQUENCE CAN BE READ PAM PROXIMAL TO PAM DISTAL
 */

unsigned long gRNA::compressSeq(std::string s) {
    unsigned long compseq;
    compseq = 0;
    for(int i=0; i<s.size(); i++) {
        unsigned long val = convertCharBase4(s[i])*pow(4,i); //multiplying by power-4 converts to base10
        compseq += val;
    }
    return compseq; //base10 version of sequence string
}


/* Function: getVectorPair
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Will take the stored sequence base-10 number and convert it to base-64 so that it may be written in the
 * output file in the smallest possible form.  This code returns a string with the total sequence compressed
 * into base-64 by combing the tailseq from the gRNA object and getting the seed sequence through a pass into
 * the function.
 */

std::pair<unsigned long, std::string> gRNA::getVectorPair(unsigned long seed) {
    //run base-10 to base-64 conversion on the location and the Sequence and then put them into loc and seq
    std::string seq = baseConvert(seed, 64);
    std::string seq2 = baseConvert(tailSeq, 64);
    std::string totseq;
    //Get absolute value of the location and add the sense signal:
    if (PAMlocation < 0) {
        PAMlocation *= -1;  //converts the location to absolute value for compression
        totseq = seq + "-" + seq2;
    } else { totseq = seq + "+" + seq2;}
    //This next part adds the score to the sequence:
    std::string scr = "," + baseConvert(OnScore,64);
    totseq += scr;
    return std::make_pair(PAMlocation, totseq);
}

/* Function: getHypLoc, getHypTail
 * --------------------------------------------------------------------------------------------------------
 * Usage: Will return a base64 string of the location with the sign designating the strand direction.
 */

std::string gRNA::getHypLoc() {
    unsigned long pam = abs(PAMlocation);
    return baseConvert(pam, 64);
}

std::string gRNA::getHypTail() {
    return baseConvert(tailSeq, 64);
}

/* Function: convertCharBase4
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Simple switch function that changes a nucleotide back and forth to base-4 representation
 */
 
int gRNA::convertCharBase4(char c) {
    switch (c) {
         case 'A': return 0;
         case 'T': return 1;
         case 'C': return 2;
         case 'G': return 3;
         default: return 0;
     }
 }

/* PUBLIC Function: baseConvert(unsigned long long base10 number, int base output)
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Simple switch from a base10 to the base specified (will be called 64 in this class).
 */

std::string gRNA::baseConvert(unsigned long long number, int base) {
    std::string base64set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789=/"; //base64 modification so +- can be used for strand direction
    std::string base32set = "ABCDEFGHIJKLMNOPQRSTUVWXYZ234567"; //base32
    int rem;
    std::string newNum;
    while(number>=base) {
        rem = number%base;
        number = number/base;
        newNum = base64set[rem] + newNum;  //only for base-64. Change the char array/string if you want to do something else.
    }
    newNum = base64set[number] + newNum;
    return newNum; //this returns the string in reverse...
}


