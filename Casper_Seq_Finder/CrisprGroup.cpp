//
//
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include <string>
#include <vector>
#include <assert.h>
#include "CrisprGroup.h"
#include "Scoring.h"
#include "gRNA.h"
#include "set.h"
#include "map.h"
#include "hashset.h"
#include "hashmap.h"
#include "foreach.h"


/* Constructor: CrisprGroup
 * Usage: CrisprGroup crisprgroup;
 * -------------------------------------------------------------------------------------------------------
 * Creates a CrisprGroup object and initializes the variables to the appropriate values
 */

CrisprGroup:: CrisprGroup(int num, std::string base, std::string org) {
    sCur = NULL;
    numChromosomes = num;
    filename = base + "NAG_files/" + org + ".txt";  // why does this say NAG_files????
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

void CrisprGroup::findPAMs (std::string &s, bool dir, int chrm, std::string pam, bool on, bool anti,std::string score_file) {  //NEED TO INCLUDE SEQUENCE LENGTH (20-24).
    //Scoring algorithm initialization:
    Scoring scoring(score_file);
    int pamsize = pam.length();
    long i;
    long index=0;
    pamEval p = evaluate(pam, anti);
    pam = p.core;
    std::cout << pam << endl;
    int offset = p.offset;
    std::vector<int> options = p.options;
    int t=0;
    while (index < s.length()) {
        i = s.find(pam,index);
        for (int k=0; k<options.size(); k++) {
            char c = s.at(i+options[k]);
            if(p.checkOptions(c,p.optid[k]) == false) {
                t++;
            }
        }
        if (i == std::string::npos) { break;}
        else if (t>0) {
            t--;
        }
        else if (s.length()-10 > i && i > 35) {
            gRNA* sequence = new gRNA;
            std::string fullseq;
            unsigned long seed;
            long j = i;
            if(!dir) { // Reverse strand detection
                j = s.length()-(i+1);
            }
            int score;
            if (anti) {
                fullseq = s.substr(i+pamsize,20);
                //scoring the full sequence including PAM and flanking regions
                score = scoring.calcScore(s.substr(i-4,pamsize+20));
                if (fullseq.length() < 20) {
                    break;
                }
                seed = sequence->insertSequence(j,chrm,pamsize,anti,dir,fullseq,score);
            } else {
                fullseq = s.substr(i-20-offset,20);
                //scoring the full sequence including PAM and flanking regions
                score = scoring.calcScore(s.substr(i-20-offset,23));
                seed = sequence->insertSequence(j,chrm,pamsize,anti,dir,fullseq,score);
            }
            //std::cout << fullseq << "," << j << std::endl; //For double checking the sequence and location
            addToMap(seed,sequence);
        }
        index = i+1;
        //reporting to determine speed of access:
        if (index%10000 == 0) {
            std::cout << index << " Positions searched." << endl;
        }
    }
}


/* Function: addToMap
 * --------------------------------------------------------------------------------------------------------
 * Usage: Checks to see if the base10 unsigned long passed into the function has already been seen by
 * comparing against the Seed_Map where all the seeds that have been seen before are stored. This
 */

void CrisprGroup::addToMap(unsigned long seed, gRNA* obj) {
    //checking if the seed sequence has already been seen
    std::unordered_map<unsigned long, std::vector<gRNA*>>::const_iterator got = Seed_Map.find(seed);
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

/* Function: evaluate
 * ------------------------------------------------------------------------------------------------------
 * Usage: Determine the PAM structure and whether the seed sequence is 5' or 3' from the PAM.
 */

pamEval CrisprGroup::evaluate(std::string pam, bool anti) {
    pamEval ret;
    if (anti) {
        std::string newpam;
        for (int i=pam.length()-1;i>=0;i--) {
            newpam += pam.at(i);
        }
        pam = newpam;
    }
    bool corefind = false;
    int pamloc=0;
    for (int i=0; i<pam.size(); i++) {
        if (pam[i] == 'N') {
            ret.offset++;
        }
        if (pam[i] == 'G' || pam[i] == 'A' || pam[i] == 'T' || pam[i] == 'C') {
            if (!corefind) {
                pamloc = i;
                corefind = true;
            }
            ret.core += pam[i];
        }
        else {
            if (pam[i] != 'N') {
                ret.options.push_back(i-pamloc);
                ret.optid.push_back(pam[i]);
            }
        }
    }
    return ret;
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
            std::pair<unsigned long, std::vector<gRNA*>> insert = std::make_pair(i.first, i.second);
            repeat_seqs.push_back(insert);
        } else { //If the sequence is unique
            gRNA* myTarget = i.second.at(0);
            //string will contain whether the location is on the sense or antisense strand
            std::pair<unsigned long, std::string> insert = myTarget->getVectorPair(i.first);
            //std::pair<long, std::string> insert = std::make_pair(myTarget->getLocation(), totalcompressed);
            total_seqs[myTarget->chrNumber()-1].push_back(insert);
            delete myTarget;
        }
    }
    for (int i=0; i<numChromosomes; i++) {
        std::cout << "sorting " << "scaffold or chromosome: " << i << "..." << std::endl;
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

















