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
    fiveSeq = 0;
    threeSeq = 0;
    //does there need to be any extra here i.e. deleting pointers?
}


/* Function: insertSequence
 * ---------------------------------------------------------------------------------------------------------
 * Usage: After an instance of gRNA is created, this function assigns the crucial parameters to the internal
 * assignments, being the complete sequence (g_sequence) and the location of the sequence stored as the
 * location of the protospacer adjacent motif (PAMlocation). This function also returns the seed sequence so
 * that it can be stored in the Seed_Map of CrisprGroup for comparison to discover repeats.
 */

unsigned long long gRNA::insertSequence(long index, int chr, int pamsize, bool sensestrand, std::string seq, int score, short five, short seedsize, short three) {
    PAMlocation = index;
	if (!sensestrand) {
        PAMlocation *= -1;
    }
    OnScore = score;
    Chromosome = chr+1; // this added 1 takes care of the indexing error and properly assigns chromosome number.
    /* The following code process the total sequence into multiple elements including the five prime and three prime sequences as well as the pam sequences. */
	if (!sensestrand) {
        pamSeq = compressSeq(seq.substr(five+seedsize+three,pamsize));
        threeSeq = compressSeq(seq.substr(five+seedsize,three));
        fiveSeq = compressSeq(seq.substr(0,five));
    } else {
        pamSeq = compressSeq(seq.substr(0,pamsize));
        fiveSeq = compressSeq(seq.substr(pamsize,five));
        threeSeq = compressSeq(seq.substr(pamsize+five+seedsize,three));
    }
    return compressSeq(seq.substr(five,seedsize));
}

/* Function: compressSeq
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Called from the insertSequence function to take in the passed sequence and compresses it into a
 * base-10 integer that can be stored with data requirements.
 * Special note: THIS PROCESS FLIPS THE STRING DIRECTION SO A SEQUENCE CAN BE READ PAM PROXIMAL TO PAM DISTAL
 */

unsigned long long gRNA::compressSeq(std::string s) {
    unsigned long compseq;
    compseq = 0;
    for(int i=0; i<s.size(); i++) {
        unsigned long val = convertCharBase4(s[i])*pow(4,i); //multiplying by power-4 converts to base10
        compseq += val;
    }
    return compseq; //base10 version of sequence string
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





