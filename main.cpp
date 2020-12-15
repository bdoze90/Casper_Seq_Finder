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
#include "WriteFile.h"
using namespace std;

/* Function Prototypes */

void selectOrganism(vector<string> &files, string &endo); //runs through the code that allows one to populate the file locations vector for analysis
string reverseComplement(string &str); //returns reverse complement of input
string toCapitals(string &str); //takes the string to all capitals


/* Main */

//the line limit for the file and the capitals mixed
//int argc, const char * argv[] -> add when exporting executable
//int main(int argc, const char * argv[]) {
int main()
{
    //int argc = 10;
	std::vector<std::string> argv = {"Executable","saCas9","NNGRRT","scede","FALSE","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERinfo","C:/Users/Tfry/Desktop/sce.fna", "Saccharomyces Cerevisiae S288C", "20", "16","notes_go_here"};
     string pamname = argv[1];
     string pam = argv[2];
     string OrgCode = argv[3];
     string returnPath = argv[5];
     bool anti = false;
     string a = string(argv[4]);
     if (a == "TRUE") {
     anti = true;
     }
    string genome_name = string(argv[8]);
    int clen = std::stoi(string(argv[9]));
    int slen = std::stoi(string(argv[10]));
    string misc_notes = string(argv[11]);
    //end obtaining information from argv.
    std::clock_t start;
    double duration;
    start = std::clock();
    string output_file = OrgCode + "_" + pamname;
    Read read;
    read.setFileName(argv[7]);
    std::cout << "Opening fasta-type file: " << argv[7] << std::endl;
    read.openFile();
    std::string score_file = argv[6];
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
    std::cout << "Finished reading in the genome file.\n";
    //Container for holding the statistics of the fasta for the end:
    std::vector<int> karystats;
    //fixes the off by one of the input sequences:
    inputSequences.push_back(newseq);
    newseq.clear();
    std::cout << inputSequences.size() << endl;
    CrisprGroup *Genome = new CrisprGroup(inputSequences.size(), returnPath, OrgCode,clen,slen);
    //Beginning of the for loop that iterates through the Fasta file to find targets
    std::cout << "Processing the genome for " << pamname << " target sequences.\n";
    for (int j=0; j<inputSequences.size(); j++) {
        string chromosomeSequence = toCapitals(inputSequences.at(j));
        inputSequences.at(j).clear();
        inputSequences.at(j).shrink_to_fit();
        karystats.push_back(chromosomeSequence.size());
        Genome->findPAMs(chromosomeSequence, true, j, pam, true, anti, score_file);
        string reverseSequence;
        reverseSequence = reverseComplement(chromosomeSequence);
        chromosomeSequence.clear();
        chromosomeSequence.shrink_to_fit();
        Genome->findPAMs(reverseSequence, false, j, pam, true, anti, score_file);
        reverseSequence.clear();
        reverseSequence.shrink_to_fit();
        cout << "Chromosome " << j+1 << " complete." << endl;
    }

    Genome->processTargets();
    cout << "Finished Locating All Cas9 target sequences" << endl;
    WriteFile Output;
    Output.inputStats(karystats, misc_notes);  // Load the statistics of the size of the chromosomes/scaffolds
    Output.setFileName(returnPath + output_file + ".cspr", genome_name);
    //Reporting the statistics:
    cout << "There were " << Genome->totSize() << " unique sequences." << endl;
    cout << "There were " << Genome->repSize() << " identical repeated sequences." << endl;
    cout << "Printing to file..." << endl;
    Output.retrieveData(Genome,chromscaff);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Time Elapsed: " << duration/60 << " minutes \n";
    delete Genome;
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
