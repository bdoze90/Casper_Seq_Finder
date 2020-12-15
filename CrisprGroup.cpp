//
//
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <assert.h>
#include "CrisprGroup.h"
#include "Scoring.h"
#include "gRNA.h"


/* Constructor: CrisprGroup
 * Usage: CrisprGroup crisprgroup;
 * -------------------------------------------------------------------------------------------------------
 * Creates a CrisprGroup object and initializes the variables to the appropriate values
 */

CrisprGroup:: CrisprGroup(int num, std::string base, std::string org, int tot_length, int seed_length) {
    sCur = NULL;
    numChromosomes = num;
    filename = base + "NAG_files/" + org + ".txt";  // why does this say NAG_files????
    len_seq = tot_length;
    len_seed = seed_length;
}

/* Destructor: ~CrisprGroup
 * Usage: delete CrisprGroup crisprgroup;
 * --------------------------------------------------------------------------------------------------------
 * Goes through the CrisprGroup object and deletes all pointer references and vectors to clear up storage
 * from the program run.
 */

CrisprGroup::~CrisprGroup() {
    //erase all the gRNA objects created at the end of main. Iterate through all of them. This should already be done by processTargets.
}

/* Function: findPAMs
 * -------------------------------------------------------------------------------------------------------
 * Usage: A recursive function that goes through the entire sequence that was inputted by the user.  It
 * finds all instances of "GG" which is the signal for a PAM sequence.  It then takes that location and
 * makes a new instance of gRNA in which the sequence is placed into to fill the data of the object.
 */

void CrisprGroup::findPAMs (std::string &s, bool dir, int chrm, std::string p, bool on, bool anti,std::string score_file) {
    //Scoring algorithm initialization:
    Scoring scoring(score_file);
    int pamsize = p.length();
    // PAM sequence search initialization
    pamEval pe;
    std::regex pam (pe.regexPAM(p));
    std::smatch m;
    auto begin = std::sregex_iterator(s.begin(), s.end(), pam);
    auto end = std::sregex_iterator();
    //searching the sequence for the PAM with the iterator
    for (std::sregex_iterator i = begin; i != end; i++) {
        std::smatch match = *i;
        std::string mpam = match.str(1);
        long m_pos = match.position();
        
        if (s.length()-10 > m_pos && m_pos > 30) { //checking to make sure you are not too close to the edge of the sequence
            gRNA* sequence = new gRNA;
            std::string fullseq;
            unsigned long seed;
            long j = m_pos;
            if(!dir) { // Reverse strand detection for appropriate location indexing on output
                j = s.length()-(m_pos+1);
            }
            int score;
            if (anti) {
                fullseq = s.substr(m_pos,len_seq+pamsize);
                //scoring the full sequence including PAM and flanking regions
                score = scoring.calcScore(s.substr(m_pos-4,pamsize+len_seq));
                //Double check to make sure the sequence is not to close to the edge so as to not get a sequence
                if (fullseq.length() < len_seq) {
                    break;
                }
                seed = sequence->insertSequence(j,chrm,pamsize,anti,dir,fullseq,score,len_seed);
            } else {
                fullseq = s.substr(m_pos-len_seq,len_seq+pamsize);
                //scoring the full sequence including PAM and flanking regions
                score = scoring.calcScore(s.substr(m_pos-len_seq,len_seq));
                seed = sequence->insertSequence(j,chrm,pamsize,anti,dir,fullseq,score,len_seed);
            }
            //std::cout << fullseq << "," << j << std::endl; //For double checking the sequence and location
            addToMap(seed,sequence);
        }
        //Reporter for how much of the sequence has been searched.
        if (m_pos%10000 == 0) {
            std::cout << m_pos << " Positions searched." << endl;
        }
        
    }
}


/* Function: addToMap
 * --------------------------------------------------------------------------------------------------------
 * Usage: Checks to see if the base10 unsigned long passed into the function has already been seen by
 * comparing against the Seed_Map where all the seeds that have been seen before are stored. This
 */

void CrisprGroup::addToMap(unsigned int seed, gRNA* obj) {
    //checking if the seed sequence has already been seen
    std::unordered_map<unsigned int, std::vector<gRNA*>>::const_iterator got = Seed_Map.find(seed);
    //if not found create a new map entry with the seed and new vector for storing any gRNA objects.
    if (got == Seed_Map.end()) {
        std::vector<gRNA*> locs;
        locs.push_back(obj);
        Seed_Map.emplace(seed,locs);
    //if it already has been seen then just add it to the vector at that map key
    } else {
        Seed_Map[seed].push_back(obj);
    }
}


/* Function: charToInt
 * -------------------------------------------------------------------------------------------------------
 * Usage: Takes in a character value representing a nucleotide and turns it into a representative integer
 */

int CrisprGroup::charToInt(char c) {
    switch (c) {
        case 'A': return 0;
        case 'T': return 1;
        case 'C': return 2;
        case 'G': return 3;
        default: return 0;
    }
}

/* Function: processTargets
 * ------------------------------------------------------------------------------------------------------
 * Usage: Takes the Seed_Map and sorts the sequences and their locations by chromosome buckets and into
 * a container for the non-unique sequences to be printed separately. Function deletes the Seed_Map at
 * the end to save memory during the program's execution.
 */

void CrisprGroup::processTargets() {
    //Initiating the total_seqs vector:
    for (int i=0; i<numChromosomes; i++) {
        std::vector<std::pair<long, std::string>> newStorVec;
        total_seqs.push_back(newStorVec);
    }
    //Iterating across the Seed_Map:
    for (const auto &i : Seed_Map) {
        if (i.second.size() != 1) { //If the sequence is non-unique
            std::pair<unsigned int, std::vector<gRNA*>> insert = std::make_pair(i.first, i.second);
            repeat_seqs.push_back(insert);
        } else { //If the sequence is unique
            gRNA* myTarget = i.second.at(0);
            //string will contain whether the location is on the sense or antisense strand
            //to generate the uncompressed sequence for debugging, set to false
            std::pair<unsigned int, std::string> insert = myTarget->getVectorPair(i.first,true);
            total_seqs[myTarget->chrNumber()-1].push_back(insert);
            delete myTarget;
        }
    }
    for (int i=0; i<numChromosomes; i++) {
        std::cout << "sorting " << "scaffold or chromosome: " << i << "..." << std::endl;
        //sort in descending order so that deletion of objects does not result in reindexing of array
        std::sort(total_seqs[i].begin(), total_seqs[i].end());
        std::cout << "done sorting." << std::endl;
    }
    Seed_Map.clear();
}

/* Function: totSize
 * ---------------------------------------------------------------------------------------------------
 * Usage: Called by main and will iterate through the chromosome vectors withing totSize to get the total
 * number of sequences.
 */

unsigned long CrisprGroup::totSize() {
    unsigned long total = 0;
    for (int i=0;i<total_seqs.size();i++) {
        total += total_seqs.at(i).size();
    }
    return total;
}


/* Function: nextUnique
 * ---------------------------------------------------------------------------------------------------
 * Usage: Called by WriteFile instance to grab the next item in the total_seqs vector. This function
 * will return a compressed (base64) string that contains the location and sequence in that order.
 */

std::string CrisprGroup::nextUnique(int chr, long index) {
    std::pair<long, std::string> current = total_seqs[chr][index];
    gRNA gRNA;
    std::string loc = gRNA.baseConvert(current.first, 64);
    std::string element = loc + "," + current.second;
    return element;
}

















