//
//  Read.cpp
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <assert.h>
using namespace std;
#include "Read.h"

void Read::openFile()
{
    stream = new std::ifstream;
    stream->open(filename, std::ifstream::in);
}


std::string Read::FirstLine()
{
    std::string junk;
    std::getline(*stream, junk);
    std::cout << "First line: " << junk << std::endl;
    return junk;
}

void Read::closeFile()
{
    stream->close();
}

string Read::getLine()
{
    std::string nts;
    std::getline(*stream, nts);
    return nts;
}

bool Read::newLine()
{
    if (!stream->eof()) return true;
    return false;
}

/* Functions to read information derived from the Genome input */

string Read::getPAM()
{
    std::string line;
    std::getline(*stream, line);
    return line.substr(17,10);
}

vector<string> Read::getOPAMs()
{
    vector<string> ret;
    std::string line;
    std::getline(*stream, line);
    std::string opam_i;
    std::getline(*stream, opam_i);
    opamnum = std::stoi(opam_i);
    for (int j = 0; j < opamnum; j++) {
        std::string aline;
        std::getline(*stream, aline);
        ret.push_back(aline);
    }
    return ret;
}

string Read::getOrgCode()
{
    std::string line;
    std::getline(*stream, line);
    return line.substr(14, 25);
}

vector<string> Read::getFileLocations()
{
    vector<string> ret;
    std::string line;
    std::getline(*stream, line);
    std::string fi;
    std::getline(*stream, fi);
    int i = std::stoi(fi);
    for (int j = 0; j < i; j++) {
        std::string aline;
        std::getline(*stream, aline);
        ret.push_back(aline);
    }
    return ret;
}

bool Read::getAnti()
{
    std::string line;
    std::getline(*stream, line);
    if (line == "TRUE") {
        return true;
    }
    else return false;
}

