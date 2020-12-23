//
//  WriteFile.cpp
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include "WriteFile.h"
using namespace std;

WriteFile::WriteFile() {
    //chromosomeseqcount = {0,0,0,0,0,0};
}

WriteFile::~WriteFile() {
    outputfile.close();
}

// See the setFileName function for incorporation of this data in the output file
void WriteFile::inputStats(std::vector<int> kary, std::string misc) {
    chromosomeseqcount = kary;
    chr_stats_str = "KARYSTATS: ";
    for (int i = 0; i<kary.size(); i++) {
        chr_stats_str += to_string(kary[i]) + ",";
    }
    mystats = "MISCELLANEOUS: " + misc;
    
}

void WriteFile::setFileName(string fn, string genome_name) {
    filename = fn;
    outputfile.open(filename);
    outputfile << "GENOME: " << genome_name << "\n";
    outputfile << chr_stats_str << "\n";
    outputfile << mystats << "\n";
}

void WriteFile::retrieveData(CrisprGroup* genome,std::vector<std::string> cs, bool repeats) {
    //retrieving the unique sequences
    std::pair<long,string> current;
    for (int i=0;i<genome->chrCount();i++) {
        outputfile << cs[i] << " (" << i+1 << ")" << "\n";
        // Loop counter is in the correct direction (positive to file).
        for (int j=0; j<genome->Size(i); j++) {
            current = genome->nextUnique(i,j);
            if (current.first < 0) {
                long outnum = (chromosomeseqcount[i] + current.first)*-1;
                outputfile << outnum << ",";
            } else {
                outputfile << current.first << ",";
            }
            outputfile << current.second << "\n";
        }
    }
    outputfile << "END_OF_FILE";
    //retrieving the repeated sequences if selected
    if (repeats) {
        //Get the pam evaluation information
        PAMstat = genome->getPamEval();
        repeatfile.open(filename + "_repeats");
        repeatfile << "REPEATS" << "\n";
        std::pair<unsigned int, std::vector<gRNA*>> newSet;
        for (int j=0;j<genome->repSize();j++) {
            newSet = genome->nextRepeatSet(j);
            string seed = decompressSeq(newSet.first,PAMstat.seedsize);
            repeatfile << seed << "\n";
            for (int i=0; i<newSet.second.size(); i++) {
                inputRepeatData(newSet.second.at(i));
                long newposition = stol(position);
                if (newposition < 0) {
                    newposition = (chromosomeseqcount[j] + newposition)*-1;
                }
                repeatfile << chromosome << "," << newposition << "," << sequence << "," << score << "\t";
                delete newSet.second.at(i);
            }
            repeatfile << "\n";
        }
    }
}

void WriteFile::inputRepeatData(gRNA* g) {
    sequence = decompressSeq(g->getFiveSeq(), PAMstat.fivesize) + "-";
    sequence += decompressSeq(g->getThreeSeq(), PAMstat.threesize);
    chromosome = g->chrNumber();
    std::string pam = decompressSeq(g->getPam(),PAMstat.pam.size());
    if (PAMstat.directionality) {
        sequence = pam + "-" + sequence;
    } else {
        sequence += "-" + pam;
    }
    score = std::to_string(g->getScore());
    position = std::to_string(g->getLocation());
}

/*void WriteFile::inputUniqueData(CrisprGroup* genome) {
 outputfile << "There are " << genome->Size() << " unique sequences across the genome. \n";
 outputfile << "There are" << genome->nagsize() << " NAG sequences across the genome. \n";
 for (int i =1; i <= chromosomeseqcount.size(); i++) {
 outputfile << "There are " << chromosomeseqcount[i] << " unique sequences on Chromosome " << i << "\n";
 }
 }*/

/* Function: decompress
 * -------------------------------------------------------------------------------------------------------
 * Usage: Takes in a long long object representing a DNA sequence and turns it into a string for printing
 */
std::string WriteFile::decompressSeq(unsigned long long cseq, short exp_len) {
    // Goes through if statement because of an off by 1 error on a zero length sequence
    if (exp_len > 0) {
        std::string uncompressed;
        //do the reverse binary transition from base-10 to base-4
        while (cseq >= 4) {
            int rem = cseq%4;
            cseq = cseq/4;
            uncompressed += convertBase4toChar(rem);
        }
        uncompressed += convertBase4toChar(cseq);
        for (int i=uncompressed.size(); i<exp_len; i++) {
            uncompressed += 'A';
        }
        return uncompressed;
    }
    return "";
}


/* Function: charToInt
 * -------------------------------------------------------------------------------------------------------
 * Usage: Takes in a character value representing a nucleotide and turns it into a representative integer
 */

int WriteFile::charToInt(char c) {
    switch (c) {
        case 'A': return 0;
        case 'T': return 1;
        case 'C': return 2;
        case 'G': return 3;
        default: return 0;
    }
}

/* Function: convertBase4toChar
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Simple switch function, reverse of above.
 */

char WriteFile::convertBase4toChar(int i) {
    std::string bfour = "ATCG";
    return bfour[i];
}

