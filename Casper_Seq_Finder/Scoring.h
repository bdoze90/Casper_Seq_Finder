//
//  Scoring.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#ifndef __Casper_Seq_Finder__Scoring__
#define __Casper_Seq_Finder__Scoring__

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <assert.h>

class Scoring {
public:
    Scoring(std::string sfile);
    double calcScore(std::string sequence);
private:
    void scanScore();
    void fillScoringAlgorithm();
private:
    std::string scoring_file;
    struct iden {
        char nt1;
        char nt2;
        int position;
        double odds_score;
    };
    std::string sequence;
    std::string PAM;
    double totalScore;
    int returnScore;
    std::vector<iden> Idens;
};

#endif /* defined(__Casper_Seq_Finder__Scoring__) */
