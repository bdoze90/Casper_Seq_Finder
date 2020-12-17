//
//  WriteFile.cpp
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include "WriteFile.h"
using namespace std;

static int callback(void *NotUsed, int argc, char **argv, char **azColName)
{
	/*
	int i;
	for (i = 0; i < argc; i++) {
		printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
	}
	printf("\n");
	*/
	return 0;
}

WriteFile::WriteFile() {
    //chromosomeseqcount = {0,0,0,0,0,0};
}

WriteFile::~WriteFile() {
    outputfile.close();
}

// See the setFileName function for incorporation of this data in the output file
void WriteFile::inputStats(std::vector<int> kary, std::string misc) {
    chr_stats_str = "KARYSTATS: ";
    for (int i = 0; i<kary.size(); i++) {
        chr_stats_str += to_string(kary[i]) + ",";
    }
    mystats = "MISCELLANEOUS: " + misc;
    
}

void WriteFile::setFileName(string fn, string genome_name) {
	filename = fn;
	outputfile.open(filename, ios_base::out | ios_base::binary);

	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
	outbuf.push(boost::iostreams::gzip_compressor());
	outbuf.push(outputfile);
	//Convert streambuf to ostream
	ostream out(&outbuf);

	out << "GENOME: " << genome_name << "\n";
	out << chr_stats_str << "\n";
	out << mystats << "\n";
}

void WriteFile::retrieveData(CrisprGroup* genome,std::vector<std::string> cs) {
    //retrieving the unique sequences
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
	outbuf.push(boost::iostreams::gzip_compressor());
	outbuf.push(outputfile);
	ostream out(&outbuf);

	std::string current;
	for (int i = 0; i < genome->chrCount(); i++) {
		out << cs[i] << " (" << i + 1 << ")" << "\n";
		// Loop counter is in the correct direction (positive to file).
		for (int j = 0; j < genome->Size(i); j++) \
		{
			current = genome->nextUnique(i, j);
			out << current << "\n";
		}
	}
	boost::iostreams::close(outbuf);
 
	//repeats
	sqlite3 *db;
	char *zErrMsg = 0;
	int rc;
	string sql;
	filename.erase(filename.find("."), 5);
	filename = filename + "_repeats.db";
	rc = sqlite3_open(filename.c_str(), &db);

	sql = "DROP TABLE IF EXISTS repeats";
	rc = sqlite3_exec(db, sql.c_str(), callback, 0, &zErrMsg);

	sql = "CREATE TABLE repeats (seed TEXT, chromosome TEXT, location TEXT, tail TEXT, score TEXT, count int);";
	rc = sqlite3_exec(db, sql.c_str(), callback, 0, &zErrMsg);

	std::pair<unsigned int, std::vector<gRNA*>> newSet;
	for (int j = 0; j < genome->repSize(); j++) 
	{
		newSet = genome->nextRepeatSet(j);
		int size = newSet.second.size();
		string seed = genome->decompressSeq(newSet.first, genome->len_seed);
		string chroms = "";
		string locs = "";
		string scores = "";
		string pams = "";
		string pam = "";
		string tail = "";
		string tails = "";
		string strand = "";

		for (int i = 0; i < newSet.second.size(); i++)
		{
			chromosome = newSet.second.at(i)->chrNumber();
			score = to_string(genome->decompress_ontarget(newSet.second.at(i)->getScore()));
			position = to_string(newSet.second.at(i)->getLocation());
			pam = genome->decompress64(newSet.second.at(i)->getHypPam(), genome->pam_length);
			tail = genome->decompress64(newSet.second.at(i)->getHypTail(), genome->len_seq - genome->len_seed);
			
			if (stol(position) >= 0)
			{
				strand = "+";
			}
			else
			{
				strand = "-";
			}

			if (i != size - 1)
			{
				chroms += to_string(chromosome) + ",";
				locs += to_string(abs(stol(position))) + ",";
				if (tail != "" && pam != "")
				{
					tails += tail + strand + pam + ",";
				}
				scores += score + ",";
			}
			else
			{
				chroms += to_string(chromosome);
				locs += to_string(abs(stol(position)));
				if (tail != "" && pam != "")
				{
					tails += tail + strand + pam;
				}
				scores += score;
			}
			delete newSet.second.at(i);
		}


		sql = "INSERT INTO repeats ('seed', 'chromosome', 'location', 'tail', 'score', 'count') VALUES (";
		sql += "'" + seed + "'" + ",";
		sql += "'" + chroms + "'" + ",";
		sql += "'" + locs + "'" + ",";
		sql += "'" + tails + "'" + ",";
		sql += "'" + scores + "'" + ",";
		sql += "'" + to_string(size) + "'";
		sql += ");";
		rc = sqlite3_exec(db, sql.c_str(), callback, 0, &zErrMsg);

	}
	sqlite3_close(db);
}



void WriteFile::inputData(gRNA* g) {
    sequence = g->getHypTail();
    chromosome = g->chrNumber();
    std::string pam = g->getHypPam();
    if (g->getLocation() < 0) {
        sequence += "-" + pam;
    } else {
        sequence += "+" + pam;
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

