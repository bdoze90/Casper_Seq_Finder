//
//  main.cpp
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include <iostream>
#include <ctime>
#include "Read.h"
#include <string.h>
#include <vector>
#include <set>
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
int main() {
    int argc = 12;
    std::vector<std::string> argv = {"none","spCas9","NGG","None","eco","FALSE","/Users/brianmendoza/Desktop/","/Users/brianmendoza/Desktop/CASPER-master/CRISPRscan.txt","/Users/brianmendoza/Dropbox/eco.fna","8","12"};
//int main(int argc, char *argv[]) {
    // argv contains in order: pamname, PAM, OPAM, OrgCode, anti, returnPath, *file locations.
     string pamname = argv[1];
     string pam = argv[2];
     vector<string> Opams;
     bool offexist = true;
     string opam = string(argv[3]);
     if (opam == "None") {
     offexist = false;
     }
     Opams.push_back(opam);
     string OrgCode = argv[4];
     string returnPath = argv[6];
     vector<string> file_locations;
     for (int l=7; l<argc; l++) {
     file_locations.push_back(argv[l]);
     }
     bool anti = false;
     string a = string(argv[5]);
     if (a == "TRUE") {
     anti = true;
     }
    int seed = std::stoi(string(argv[9]));
    int tail = std::stoi(string(argv[10]));
    //end obtaining information from argv.
    std::clock_t start;
    double duration;
    start = std::clock();
    string output_file = OrgCode + pamname;
    Read read;
    read.setFileName(argv[8]);
    std::cout << "Opening fasta-type file: " << argv[8] << std::endl;
    read.openFile();
    std::string score_file = argv[7];
    //input sequences need to be a vector...
    vector<string> inputSequences;
    string newseq = "";
    read.skipLine(); //skips the first line of the title of the fasta file
    std::cout << "Reading file...\n";
    while (read.newLine()) {
        std::string line = read.getLine(); //WARNING: THIS ONLY ACCOMODATES A 100 NUCLEOTIDE LINE!
        if (line[0] == '>') {
            std::cout << "Chromosome type " << line << " detected.\n";
            inputSequences.push_back(newseq);
            newseq = "";
            while (line[0] != '>') {
                line = read.getLine();
            }
        } else {
            newseq += line; //THIS ACCOMODATES UP TO A 100000000 NUCLEOTIDE LINE
        }
    }
    std::cout << "Finished reading in the genome file.\n";
    //fixes the off by one of the input sequences:
    inputSequences.push_back(newseq);
    newseq.clear();
    std::cout << inputSequences.size() << endl;
    CrisprGroup Genome(inputSequences.size(), returnPath, OrgCode);
    //Beginning of the for loop that iterates through the Fasta file to find targets
    std::cout << "Processing the genome for " << pamname << " target sequences.\n";
    for (int j=0; j<inputSequences.size(); j++) {
        string chromosomeSequence = toCapitals(inputSequences.at(j));
        inputSequences.at(j).clear();
        inputSequences.at(j).shrink_to_fit();
        Genome.findPAMs(chromosomeSequence, true, j, pam, true, anti, score_file,seed,tail);
        if (offexist) {
            Genome.findPAMs(chromosomeSequence, true, j, Opams[0], false, anti, score_file,seed,tail);
        }
        string reverseSequence;
        reverseSequence = reverseComplement(chromosomeSequence);
        chromosomeSequence.clear();
        chromosomeSequence.shrink_to_fit();
        Genome.findPAMs(reverseSequence, false, j, pam, true, anti, score_file,seed,tail);
        if (offexist) {
            Genome.findPAMs(reverseSequence, false, j, Opams[0], false, anti, score_file,seed,tail);
        }
        reverseSequence.clear();
        reverseSequence.shrink_to_fit();
        cout << "Chromosome " << j+1 << " complete." << endl;
    }

    Genome.processTargets();
    cout << "Finished Locating All Cas9 target sequences" << endl;
    cout << "Printing to file..." << endl;
    WriteFile Output;
    Output.setFileName(returnPath + output_file + ".cspr");
    Output.retrieveData(&Genome);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Time Elapsed: " << duration << "\n";
    //Reporting the statistics:
    cout << "There were " << Genome.totSize() << " unique sequences." << endl;
    cout << "There were " << Genome.repSize() << " identical repeated sequences." << endl;
    cout << "Finished Creating File.\n To search restart CASPER and select Organism." << endl;
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
