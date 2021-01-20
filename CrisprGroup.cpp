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
#include "RepStor.h"

/* Constructor: CrisprGroup
 * Usage: CrisprGroup crisprgroup;
 * -------------------------------------------------------------------------------------------------------
 * Creates a CrisprGroup object and initializes the variables to the appropriate values
 */

CrisprGroup:: CrisprGroup(unsigned long num, pamEval mypam, std::string base, std::string org, bool repeats) {
    sCur = NULL;
    numChromosomes = num;
    filename = base + "NAG_files/" + org + ".txt";  // why does this say NAG_files????
    PAMstat = mypam;
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

/* Function: initiateTotalSeqs
* --------------------------------------------------------------------------------------------------------
* Initiate the total seqs vector so that it has allocated the appropriate vector memory.
*/

void CrisprGroup::initiateTotalSeqs() {
    //Initiating the total_seqs vector:
    for (int i=0; i<numChromosomes; i++) {
        std::vector<std::pair<long, compgrna>> newStorVec;
        total_seqs.push_back(newStorVec);
    }
}

/* Function: findPAMs
 * -------------------------------------------------------------------------------------------------------
 * Usage: A function that goes through the entire sequence that was inputted by the user.  It
 * finds all instances of "GG" which is the signal for a PAM sequence.  It then takes that location and
 * makes a new instance of gRNA in which the sequence is placed into to fill the data of the object.
 */

void CrisprGroup::findPAMs (std::string &s, bool strand, int chrm, std::string score_file) {
    //Scoring algorithm initialization:
    Scoring scoring(score_file);
    // PAM sequence search initialization
    short pamsize = PAMstat.pam.length();
    short fulllen = PAMstat.fulllen();
    std::regex pam (PAMstat.regexPAM());
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
            if(!strand) { // Reverse strand detection for appropriate location indexing on output
                j = s.length()-(m_pos+1);
            }
            if (PAMstat.directionality) {
                fullseq = s.substr(m_pos,fulllen+pamsize);
                //Double check to make sure the sequence is not to close to the edge so as to not get a sequence
                if (fullseq.length() < PAMstat.fulllen()) {
                    break;
                }
            } else {
                fullseq = s.substr(m_pos-fulllen,fulllen+pamsize);
            }
            seed = sequence->insertSequence(j,chrm,pamsize,strand,PAMstat.directionality,fullseq,scoring.calcScore(fullseq),PAMstat.fivesize,PAMstat.seedsize,PAMstat.threesize);
            //std::cout << fullseq << "," << j << std::endl; //For double checking the sequence and location
            addToMap(seed,sequence,true);
        }
        //Reporter for how much of the sequence has been searched.
        /*if (m_pos%10000 == 0) {
            std::cout << m_pos << " Positions searched." << endl;
        }*/
        
    }
}

/* Function: findPAMs_notRepeats
 * -------------------------------------------------------------------------------------------------------
 * Usage: Function that goes through the entire sequence and finds instances of the PAM sequence inputted.
 * Stores directly into the total_seqs vector.
 */

void CrisprGroup::findPAMs_notRepeats(std::string &s, bool strand, int chrm, std::string score_file) {
    //Scoring algorithm initialization:
    Scoring scoring(score_file);
    pamEval pe;
    short pamsize = pe.pam.length();
    short fulllen = pe.fulllen();
    std::regex pam (pe.regexPAM());
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
            if(!strand) { // Reverse strand detection for appropriate location indexing on output
                j = s.length()-(m_pos+1);
            }
            if (pe.directionality) {
                fullseq = s.substr(m_pos,fulllen+pamsize);
                //Double check to make sure the sequence is not to close to the edge so as to not get a sequence
                if (fullseq.length() < pe.fulllen()) {
                    break;
                }
            } else {
                fullseq = s.substr(m_pos-fulllen,fulllen+pamsize);
            }
            seed = sequence->insertSequence(j,chrm,pamsize,strand,PAMstat.directionality,fullseq,scoring.calcScore(fullseq),pe.fivesize,pe.seedsize,pe.threesize);
            //std::cout << fullseq << "," << j << std::endl; //For double checking the sequence and location
            // Transfer objects straight to the output vector to avoid sorting:
            addToMap(seed,sequence,false);
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
 * comparing against the Seed_Map where all the seeds that have been seen before are stored.  This function
 * will only add gRNA objects to the "unique seeds" vector if the user has selected not to analyze for repeats.
 */

void CrisprGroup::addToMap(unsigned long seed, gRNA* obj, bool repeats) {
    if (repeats) {
        //checking if the seed sequence has already been seen
        std::unordered_map<unsigned long, std::vector<gRNA*>>::const_iterator got = Seed_Map.find(seed);
        //if not found create a new map entry with the seed and new vector for storing any gRNA objects.
        if (got == Seed_Map.end()) {
            std::vector<gRNA*> locs;
            locs.reserve(2);
            locs.push_back(obj);
            Seed_Map.emplace(seed,locs);
            //if it already has been seen then just add it to the vector at that map key
        } else {
            Seed_Map[seed].push_back(obj);
        }
    } /*else {
        total_seqs[obj->chrNumber()].push_back();
    }*/
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
 * the end to save memory during the program's execution.  This function is not called when repeats are
 * not analyzed.
 */

void CrisprGroup::processTargets(std::vector<int> kary) {
    //Iterating across the Seed_Map:
    for (const auto &i : Seed_Map) {
        if (i.second.size() != 1) { //If the sequence is non-unique
            //unsigned int is the seed sequence and the gRNA pointer contains the rest of the info
            std::pair<unsigned long, std::vector<gRNA*>> insert = std::make_pair(i.first, i.second);
            repeat_seqs.push_back(insert);
        } else { //If the sequence is unique
            gRNA* myTarget = i.second.at(0);
			long revised_location = myTarget->getLocation();
            compgrna mystruct;
            mystruct.cfive = myTarget->getFiveSeq();
            mystruct.cthree = myTarget->getThreeSeq();
            mystruct.cpam = myTarget->getPam();
            mystruct.score = myTarget-> getScore();
            mystruct.seed = i.first;
			if (revised_location < 0) {
				revised_location = kary[myTarget->chrNumber() - 1] + revised_location;
				mystruct.strand = false;
			}
			std::pair<long, compgrna> insert = std::make_pair(revised_location, mystruct);
            //string will contain whether the location is on the sense or antisense strand
            total_seqs[myTarget->chrNumber()-1].push_back(insert);
            delete myTarget;
        }
    }
    for (int i=0; i<numChromosomes; i++) {
        std::cout << "sorting " << "scaffold or chromosome: " << i << "..." << std::endl;
        //sort in descending order so that deletion of objects does not result in reindexing of array
        std::sort(total_seqs[i].begin(), total_seqs[i].end(),pairCompare);
        //std::cout << "done sorting." << std::endl;
    }
    Seed_Map.clear();
}

bool CrisprGroup::pairCompare(const std::pair<long, compgrna> i, const std::pair<long, compgrna> j) {
    return i.first < j.first;
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

std::pair<long, std::string> CrisprGroup::nextUnique(int chr, long index) {
    std::pair <long, compgrna> cur = total_seqs[chr][index];
	std::string output_element = decompressSeq(cur.second.cfive, PAMstat.fivesize) + decompressSeq(cur.second.seed, PAMstat.seedsize) + decompressSeq(cur.second.cthree, PAMstat.threesize) + "," + decompressSeq(cur.second.cpam, PAMstat.pam.size()) + "," + std::to_string(cur.second.score);
	std::pair<long, std::string> output;
	output.first = cur.first;
	if (!cur.second.strand) {
		output.first *= -1;
	}
	output.second = output_element;
	return output;
}

/* Function: decompressSeq
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Takes the stored sequence (which is compressed to base-10) and outputs it as its string of A,T,C,G.
 * There is a problem with an Adenine in the first part of the seedSeq or tailSeq where it is not taken into
 * account because of the fact that it is represented by a 0.  To remedy this, when the decompression occurs
 * we need to make sure that if the decompression reveals a sequence shorter than the expected length, to add
 * an A to the end.
 */

std::string CrisprGroup::decompressSeq(unsigned long cseq, short exp_len) {
	// Goes through if statement because of an off by 1 error on a zero length sequence
	if (exp_len > 0) {
		std::string uncompressed;
		//do the reverse binary transition from base-10 to base-4
		while (cseq >= 4) {
			int rem = cseq % 4;
			cseq = cseq / 4;
			uncompressed += convertBase4toChar(rem);
		}
		uncompressed += convertBase4toChar(cseq);
		for (int i = uncompressed.size(); i < exp_len; i++) {
			uncompressed += 'A';
		}
		return uncompressed;
	}
	return "";
}


/* Function: convertBase4toChar
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Simple switch function, reverse of above.
 */

char CrisprGroup::convertBase4toChar(int i) {
    std::string bfour = "ATCG";
    return bfour[i];
}

/* PUBLIC Function: baseConvert(unsigned long long base10 number, int base output)
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Simple switch from a base10 to the base specified.
 */

std::string CrisprGroup::baseConvert(unsigned long long number, int base) {
    int rem;
    std::string newNum;
    while(number>=base) {
        rem = number%base;
        number = number/base;
    }
    return newNum; //this returns the string in reverse...
}

















