//
//  main.cpp
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include "Read.h"
#include <string.h>
#include <vector>
#include <set>
#include "pameval.h"
#include "gRNA.h"
#include "CrisprGroup.h"
#include "set.h"
#include "WriteFile.h"
using namespace std;

/* Function Prototypes */

void selectOrganism(vector<string> &files, string &endo); //runs through the code that allows one to populate the file locations vector for analysis
string reverseComplement(string &str); //returns reverse complement of input
string toCapitals(string &str); //takes the string to all capitals


/* Main */

//the line limit for the file and the capitals mixed
//int argc, const char * argv[] -> add when exporting executable
int main() {
    //int argc = 10;
    std::vector<std::string> argv = {"Executable","saCas9","NNGRRT", "TRUE", "FALSE", "4","16","0","testfile","/Users/brianmendoza/Desktop/","/Users/brianmendoza/Dropbox/CASPER/CASPERinfo","/Users/brianmendoza/Downloads/sceR64.fna", "Saccharomyces cerevisae", "notes_go_here"};
    pamEval P;
    P.PAMID = argv[1];
    P.pam = argv[2];
    if (string(argv[4]) == "TRUE") {
        P.directionality = true;
    } else {
        P.directionality = false;
    }
    P.fivesize = stoi(string(argv[5]));
    P.seedsize = stoi(string(argv[6]));
    P.threesize = stoi(string(argv[7]));
    
    bool repeats = false;
    string r = string(argv[3]);
    if (r == "TRUE") {
        repeats = true;
    }
    string OrgCode = argv[8];
    string returnPath = argv[9];
    string genome_name = string(argv[12]);
    string misc_notes = string(argv[13]);
    //end obtaining information from argv.
    std::clock_t start;
    double duration;
    start = std::clock();
    string output_file = OrgCode + "_" + P.PAMID;
    Read read;
    read.setFileName(argv[11]);
    std::cout << "Opening fasta-type file: " << argv[11] << std::endl;
    read.openFile();
    std::string score_file = argv[10];
    //input sequences need to be a vector...
    vector<string> inputSequences;
    string newseq = "";
    std::cout << "Reading file...\n";
    std::vector<std::string> chromscaff;
    chromscaff.push_back(read.FirstLine());  //reports the first line of the title of the fasta file and adds it to the chromscaff
    while (read.newLine()) {
        std::string line = read.getLine();
        if (line[0] == '>') {
            std::cout << "New Chromosome/Scaffold detected.\n";
            chromscaff.push_back(line);
            inputSequences.push_back(newseq);
            newseq = "";
        } else {
            newseq += line; //THIS ACCOMODATES UP TO A 100000000 NUCLEOTIDE LINE
        }
    }
    std::cout << "Finished reading in the fasta file.\n";
    //Container for holding the statistics of the fasta for the end:
    std::vector<int> karystats;
    //fixes the off by one of the input sequences:
    inputSequences.push_back(newseq);
    newseq.clear();
    std::cout << "Number of Chromosomes/Scaffolds: " << inputSequences.size() << endl;
    CrisprGroup *Genome = new CrisprGroup(inputSequences.size(), P, returnPath, OrgCode, repeats);
    //Beginning of the for loop that iterates through the Fasta file to find targets
    std::cout << "Processing the genome for " << P.PAMID << " target sequences.\n";
    Genome->initiateTotalSeqs();
    for (int j=0; j<inputSequences.size(); j++) {
        string chromosomeSequence = toCapitals(inputSequences.at(j));
        inputSequences.at(j).clear();
        inputSequences.at(j).shrink_to_fit();
        karystats.push_back(chromosomeSequence.size());
        if (repeats) {
            Genome->findPAMs(chromosomeSequence, true, j, score_file);
        } else {
            cout << "Ignoring 'repeats comparison" << endl;
            Genome->findPAMs_notRepeats(chromosomeSequence, true, j, score_file);
        }
        string reverseSequence;
        reverseSequence = reverseComplement(chromosomeSequence);
        chromosomeSequence.clear();
        chromosomeSequence.shrink_to_fit();
        if (repeats) {
            Genome->findPAMs(reverseSequence, false, j, score_file);
        } else {
            cout << "Ignoring repeats comparison" << endl;
            Genome->findPAMs_notRepeats(reverseSequence,false, j, score_file);
        }
        reverseSequence.clear();
        reverseSequence.shrink_to_fit();
        cout << "Chromosome/Scaffold " << j+1 << " complete." << endl;
    }
    if (repeats) {
        Genome->processTargets();
    }
    //std::remove(argv[8].c_str());
    cout << "Deleted temporary file" << endl;
    cout << "Finished Locating All Cas9 target sequences" << endl;
    WriteFile Output;
    Output.inputStats(karystats, misc_notes);  // Load the statistics of the size of the chromosomes/scaffolds
    Output.setFileName(returnPath + output_file + ".cspr", genome_name);
    //Reporting the statistics:
    cout << "There were " << Genome->totSize() << " unique sequences." << endl;
    cout << "There were " << Genome->repSize() << " identical repeated sequences." << endl;
    cout << "Printing to file..." << endl;
    Output.retrieveData(Genome,chromscaff,repeats);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Time Elapsed: " << duration/60 << " minutes \n";
    delete Genome;
    cout << "Finished Creating File.\n To search return to CASPER homepage and select " << genome_name << " to find targets." << endl;
    return 0;
}

/* Function: toCapitals
 * --------------------------------------------------------------------------------------------------------
 * Usage: Takes the input sequence and changes all instances of the lower case into an upper case to be passed
 * into the program for editing.
 */

string toCapitals(string &str) {
    string tc = "";
    for (int i=0; i<str.length(); i++) {
        char n = str.at(i);
        char change;
        switch (n) {
            case 'a': change = 'A'; break;
            case 't': change = 'T'; break;
            case 'g': change = 'G'; break;
            case 'c': change = 'C'; break;
            default: change = n;
        }
        tc += change;
    }
    return tc;
}

/* Function: reverseComplement
 * --------------------------------------------------------------------------------------------------------
 * Usage: a sequence in the form of a string is passed in by reference and the function returns the reverse
 * complement of the passed in sequence, inserting X's if there are any nucleotide discrepancies in the
 * original sequence.
 */

string reverseComplement(string &str) {
    string rc = "";
    for (long i=str.length()-1; i >= 0; i--) {
        char n = str.at(i);
        char reverse;
        switch (n) {
            case 'A': reverse = 'T'; break;
            case 'T': reverse = 'A'; break;
            case 'G': reverse = 'C'; break;
            case 'C': reverse = 'G'; break;
            default: reverse = 'N';
        }
        rc += reverse;
    }
    return rc;
}
