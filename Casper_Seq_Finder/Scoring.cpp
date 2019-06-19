//
//  Scoring.cpp
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include "Scoring.h"

/* Constructor for Scoring class object. One object used per CrisprGroup NOT USED IN SEQ FINDER USE IN PYTHON CODE!!! */

Scoring::Scoring (std::string sfile) {
    scoring_file = sfile;
    fillScoringAlgorithm();
}

double Scoring::calcScore(std::string s) {
    totalScore = 0;
    returnScore = 0;
    sequence = s;
    scanScore();
    //following line normalizes to best possible score
    totalScore = 1-((1.29401-totalScore)/1.94947);
    returnScore = (totalScore*100)+0.5; //0.5 for proper rounding
    return returnScore;
}

/* Iterates across the CRISPRscan scores and looks for the features in the sequence */
void Scoring::scanScore() {
    for (int i=0; i<Idens.size(); i++) {
        char nucleo1 = Idens.at(i).nt1;  //may have to change this to pointers b/c of multiple creations of object
        char nucleo2 = Idens.at(i).nt2;
        int pos = Idens.at(i).position;
        if (pos < sequence.size()) {
            if (nucleo2 != 'x') {
                std::string dinucleo = std::string() + nucleo1 + nucleo2;
                if (sequence.substr(pos-1,2) == dinucleo) {
                    totalScore += Idens.at(i).odds_score;
                }
            } else {
                if (sequence.at(pos-1) == nucleo1) {
                    totalScore += Idens.at(i).odds_score;
                }
            }
        }
    }
}

void Scoring::fillScoringAlgorithm() {
    //Establish file stream
    std::ifstream* stream = new std::ifstream;
    stream->open(scoring_file.c_str(), std::ifstream::in);
    
    //Get to the CRISPRSCAN data section of the file:
    std::string myline = "myline";
    while (myline != "CRISPRSCAN_DATA") {
        std::getline(*stream, myline);
    }
    
    //Load Information
    std::string nts;
    std::getline(*stream,nts);
    while (nts[0] != '-') {
        iden nid;
        std::vector<std::string> mytoke = Msplit(nts, '\t');
        nid.nt1 = mytoke[0][0];
        nid.nt2 = mytoke[0][1];
        nid.position = std::stoi(mytoke[1]);
        nid.odds_score = std::stod(mytoke[2]);
        Idens.push_back(nid);
        std::getline(*stream,nts);
    }
    stream->close();
}

/* A method for splitting the repeats section. */
std::vector<std::string> Scoring::Msplit(const std::string &text, char sep) {
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}


