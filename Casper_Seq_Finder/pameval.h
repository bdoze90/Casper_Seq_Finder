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

class pamEval {
public:
    std::string core;
    int offset = 0;
    std::vector<char> optid;
    std::vector<int> options;
    bool checkOptions(char n, char c) {
        bool check = false;
        switch (c) {
            case 'W':
                if (n == 'A' || n == 'T')
                    check = true;
                break;
            case 'S':
                if (n == 'C' || n == 'G')
                    check = true;
                break;
            case 'M':
                if (n == 'A' || n == 'C')
                    check = true;
                break;
            case 'K':
                if (n == 'G' || n == 'T')
                    check = true;
                break;
            case 'R':
                if (n == 'A' || n == 'G')
                    check = true;
                break;
            case 'Y':
                if (n == 'C' || n == 'T')
                    check = true;
                break;
            case 'B':
                if (n != 'A')
                    check = true;
                break;
            case 'D':
                if (n != 'C')
                    check = true;
                break;
            case 'H':
                if (n != 'G')
                    check = true;
                break;
            case 'V':
                if (n != 'T')
                    check = true;
                break;
        }
        return check;
    }
};


#endif
