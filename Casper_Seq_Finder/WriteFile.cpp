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

void WriteFile::setFileName(string fn) {
    filename = fn;
    outputfile.open(filename);
}

void WriteFile::retrieveData(CrisprGroup* genome) {
    //retrieving the unique sequences
    std::string current;
    for (int i=0;i<genome->chrCount();i++) {
        outputfile << "CHROMOSOME #" << i+1 << "\n";
        for (int j=0; j<genome->Size(i); j++) {
            current = genome->nextUnique(i,j);
            outputfile << current << "\n";
        }
    }
    //retrieving the repeated sequences
    outputfile << "REPEATS" << "\n";
    gRNA grna;
    std::pair<unsigned long, std::vector<gRNA*>> newSet;
    for (int j=0;j<genome->repSize();j++) {
        newSet = genome->nextRepeatSet(j);
        string seed = grna.baseConvert(newSet.first, 64);
        outputfile << seed << "\n";
        for (int i=0; i<newSet.second.size(); i++) {
            inputData(newSet.second.at(i));
            outputfile << chromosome << "," << position << "," << sequence << "," << score << "\t";
            delete newSet.second.at(i);
        }
        outputfile << "\n";
    }
    outputfile << "END_OF_FILE";
}

void WriteFile::inputData(gRNA* g) {
    sequence = g->getHypTail();
    chromosome = g->chrNumber();
    if (g->getLocation() < 0) {
        sequence = "-" + sequence;
    } else {
        sequence = "+" + sequence;
    }
    score = g->getScore();
    position = g->getHypLoc();
}

/*void WriteFile::printInfo(CrisprGroup* genome) {
 outputfile << "There are " << genome->Size() << " unique sequences across the genome. \n";
 outputfile << "There are" << genome->nagsize() << " NAG sequences across the genome. \n";
 for (int i =1; i <= chromosomeseqcount.size(); i++) {
 outputfile << "There are " << chromosomeseqcount[i] << " unique sequences on Chromosome " << i << "\n";
 }
 
 
 }*/

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

