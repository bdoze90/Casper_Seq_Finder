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
    stream = fopen( filename.c_str(), "r" );
    assert(stream);
}

double Read::getDouble()
{
    double rv;
    fscanf( stream, "%lf\n", &rv );
    return rv;
}

int Read::getInt()
{
    int rv;
    fscanf( stream, "%d\n", &rv);
    return rv;
}


void Read::skipLine()
{
    char junk[150];
    fscanf( stream, "%149[^\n]", junk );
    std::cout << "First line: " << string(junk) << std::endl;
}

void Read::closeFile()
{
    if(stream) fclose( stream );
}

string Read::getLine()
{
    char* nts = new char[100000000];
    fscanf( stream, "%s\n", nts );
    string retstr = string(nts);
    delete[] nts;
    return retstr;
}

bool Read::newLine()
{
    if (feof(stream) == 0 ) return true;
    return false;
}

/* Functions to read information derived from the Genome input */

string Read::getPAM()
{
    char line[25];
    fscanf(stream, "%s\n", line );
    string retstr = string(line).substr(17,10);
    return retstr;
}

vector<string> Read::getOPAMs()
{
    vector<string> ret;
    char line[20];
    fscanf(stream, "%s\n", line );
    int i;
    fscanf(stream, "%d\n", &i);
    opamnum = i;
    for (int j=0; j<i; j++) {
        char aline[10];
        fscanf(stream, "%s\n", aline);
        ret.push_back(string(aline));
    }
    return ret;
}

string Read::getOrgCode()
{
    char line[25];
    fscanf(stream, "%s\n", line);
    return string(line).substr(14,25);
}

vector<string> Read::getFileLocations()
{
    vector<string> ret;
    char line[20];
    fscanf(stream, "%s\n", line );
    int i;
    fscanf(stream, "%d\n", &i);
    for (int j=0; j<i; j++) {
        char aline[100];
        fscanf(stream, "%s\n", aline);
        ret.push_back(string(aline));
    }
    return ret;
}

bool Read::getAnti()
{
    char line[10];
    fscanf(stream, "%5s", line);
    string anti = string(line);
    if (anti == "TRUE") {
        return true;
    } else return false;
}

