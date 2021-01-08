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
	chromosomeseqcount = kary;
	chr_stats_str = "KARYSTATS: ";
    for (int i = 0; i<kary.size(); i++) {
        chr_stats_str += to_string(kary[i]) + ",";
    }
    mystats = "MISCELLANEOUS: " + misc;
    
}

void WriteFile::setFileName(string fn, string genome_name) {
	filename = fn;
	outputfile.open(filename, ios_base::out | ios_base::binary);

	//ofstream outputfile("test.cspr", ios_base::out | ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
	outbuf.push(boost::iostreams::gzip_compressor());
	outbuf.push(outputfile);
	//Convert streambuf to ostream
	ostream out(&outbuf);

	out << "GENOME: " << genome_name << "\n";
	out << chr_stats_str << "\n";
	out << mystats << "\n";
}

void WriteFile::retrieveData(CrisprGroup* genome, std::vector<std::string> cs, bool repeats) {
	//retrieving the unique sequences
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
	outbuf.push(boost::iostreams::gzip_compressor());
	outbuf.push(outputfile);
	//Convert streambuf to ostream
	ostream out(&outbuf);
	std::pair<long, string> current;
	for (int i = 0; i < genome->chrCount(); i++) 
	{
		out << cs[i] << " (" << i + 1 << ")" << "\n";
		// Loop counter is in the correct direction (positive to file).
		for (int j = 0; j < genome->Size(i); j++) 
		{
			current = genome->nextUnique(i, j);
			out << current.first << "," << current.second << "\n";
		}
	}

	boost::iostreams::close(outbuf);
	//retrieving the repeated sequences if selected

	if (repeats) {
		//setup db file
		sqlite3 *db;
		char *zErrMsg = 0;
		int rc;
		string sql;
		filename.erase(filename.find("."), 5);
		filename = filename + "_repeats.db";
		rc = sqlite3_open(filename.c_str(), &db);

		sql = "DROP TABLE IF EXISTS repeats";
		rc = sqlite3_exec(db, sql.c_str(), callback, 0, &zErrMsg);

		sql = "CREATE TABLE repeats (seed TEXT, chromosome TEXT, location TEXT, three TEXT, five TEXT, pam TEXT, score TEXT, count int);";
		rc = sqlite3_exec(db, sql.c_str(), callback, 0, &zErrMsg);

		rc = sqlite3_exec(db, "BEGIN TRANSACTION;", callback, 0, &zErrMsg);

		//Get the pam evaluation information
		PAMstat = genome->getPamEval();

		std::pair<unsigned int, std::vector<gRNA*>> newSet;
		for (int j = 0; j < genome->repSize(); j++) 
		{
			newSet = genome->nextRepeatSet(j);
			int size = newSet.second.size();
			string seed = decompressSeq(newSet.first, PAMstat.seedsize);
			string chroms = "";
			string locs = "";
			string scores = "";
			string threes = "";
			string fives = "";
			string pams = "";
			string three = "";
			string five = "";
			string pam = "";

			for (int i = 0; i < newSet.second.size(); i++)
			{
				chromosome = newSet.second.at(i)->chrNumber();
				score = to_string(newSet.second.at(i)->getScore());
				position = to_string(newSet.second.at(i)->getLocation());
				three = decompressSeq(newSet.second.at(i)->getThreeSeq(), PAMstat.threesize);
				five = decompressSeq(newSet.second.at(i)->getFiveSeq(), PAMstat.fivesize);
				pam = decompressSeq(newSet.second.at(i)->getPam(), PAMstat.pam.size());
				long newposition = stol(position);
				if (newposition < 0) 
				{
					newposition = (chromosomeseqcount[chromosome - 1] + newposition)*-1;
				}
				position = to_string(newposition);
				//inputRepeatData(newSet.second.at(i));
				//repeatfile << chromosome << "," << position << "," << sequence << "," << score << "\t";
				if (i != size - 1)
				{
					chroms += to_string(chromosome) + ",";
					locs += position + ",";
					if (three != "")
					{
						threes += three + ",";
					}
					if (five != "")
					{
						fives += five + ",";
					}
					if (pam != "")
					{
						pams += pam + ",";
					}
					scores += score + ",";

				}
				else
				{
					chroms += to_string(chromosome);
					locs += position;
					if (three != "")
					{
						threes += three;
					}
					if (five != "")
					{
						fives += five;
					}
					if (pam != "")
					{
						pams += pam;
					}
					scores += score;
				}
				delete newSet.second.at(i);
			}



			sql = "INSERT INTO repeats ('seed', 'chromosome', 'location', 'three', 'five', 'pam', 'score', 'count') VALUES (";
			sql += "'" + seed + "'" + ",";
			sql += "'" + chroms + "'" + ",";
			sql += "'" + locs + "'" + ",";
			sql += "'" + threes + "'" + ",";
			sql += "'" + fives + "'" + ",";
			sql += "'" + pams + "'" + ",";
			sql += "'" + scores + "'" + ",";
			sql += "'" + to_string(size) + "'";
			sql += ");";
			rc = sqlite3_exec(db, sql.c_str(), callback, 0, &zErrMsg);

		}

		rc = sqlite3_exec(db, "END TRANSACTION;", callback, 0, &zErrMsg);

		sqlite3_close(db);
	}

}

void WriteFile::inputRepeatData(gRNA* g) {
    sequence = decompressSeq(g->getFiveSeq(), PAMstat.fivesize) + "-";
    sequence += decompressSeq(g->getThreeSeq(), PAMstat.threesize);
    chromosome = g->chrNumber();
    std::string pam = decompressSeq(g->getPam(),PAMstat.pam.size());
    if (PAMstat.directionality) {
        sequence = pam + "-" + sequence;
    } else {
        sequence += "-" + pam;
    }
    score = std::to_string(g->getScore());
    position = std::to_string(g->getLocation());
}

/*void WriteFile::inputUniqueData(CrisprGroup* genome) {
 outputfile << "There are " << genome->Size() << " unique sequences across the genome. \n";
 outputfile << "There are" << genome->nagsize() << " NAG sequences across the genome. \n";
 for (int i =1; i <= chromosomeseqcount.size(); i++) {
 outputfile << "There are " << chromosomeseqcount[i] << " unique sequences on Chromosome " << i << "\n";
 }
 }*/

/* Function: decompress
 * -------------------------------------------------------------------------------------------------------
 * Usage: Takes in a long long object representing a DNA sequence and turns it into a string for printing
 */
std::string WriteFile::decompressSeq(unsigned long long cseq, short exp_len) {
    // Goes through if statement because of an off by 1 error on a zero length sequence
    if (exp_len > 0) {
        std::string uncompressed;
        //do the reverse binary transition from base-10 to base-4
        while (cseq >= 4) {
            int rem = cseq%4;
            cseq = cseq/4;
            uncompressed += convertBase4toChar(rem);
        }
        uncompressed += convertBase4toChar(cseq);
        for (int i=uncompressed.size(); i<exp_len; i++) {
            uncompressed += 'A';
        }
        return uncompressed;
    }
    return "";
}


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

/* Function: convertBase4toChar
 * ---------------------------------------------------------------------------------------------------------
 * Usage: Simple switch function, reverse of above.
 */

char WriteFile::convertBase4toChar(int i) {
    std::string bfour = "ATCG";
    return bfour[i];
}

