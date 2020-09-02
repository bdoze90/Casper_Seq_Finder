//
//  pameval.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#ifndef Casper_Seq_Finder_pameval_h
#define Casper_Seq_Finder_pameval_h

#include <vector>
#include <string>

class pamEval {
public:
    std::string PAMID;
    std::string pam;
    short head;
    short tail;
    short seed;
    bool directionality;
public:
    short fulllen () {
        return head + seed + tail;
    };
    std::string regexPAM() {
        std::string retpam = "(?=(";  //initializes lookahead structure
        // Iterate across the pam and generate regex characters for degenerates
        for(int i=0;i<pam.size();i++) {
            if (pam[i] == 'N') {
                retpam += '.';
            } else if (pam[i] == 'A' || pam[i] == 'T' || pam[i] == 'C' || pam[i] == 'G') {
                retpam += pam[i];
            } else {
                retpam += degenerateRegex(pam[i]);
            }
        }
        // finish the lookahead regex structure:
        retpam += ")).";
        return retpam;
    };
private:
    // Function for assigning a regex code for a degenerate nucleotide code
    std::string degenerateRegex(char c) {
        std::string reg_s;
        switch (c) {
            case 'W':
                reg_s = "[AT]";
                break;
            case 'S':
                reg_s = "[CG]";
                break;
            case 'M':
                reg_s = "[AC]";
                break;
            case 'K':
                reg_s = "[GT]";
                break;
            case 'R':
                reg_s = "[AG]";
                break;
            case 'Y':
                reg_s = "[CT]";
                break;
            case 'B':
                reg_s = "[TCG]";
                break;
            case 'D':
                reg_s = "[ATG]";
                break;
            case 'H':
                reg_s = "[ATC]";
                break;
            case 'V':
                reg_s = "[AGC]";
                break;
        }
        return reg_s;
    };
};


#endif
