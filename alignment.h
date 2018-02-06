/*
* Michelle Farr
*/
#include <iostream>
#include <string>
#include <fstream>

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::fstream;

struct DP_cell{
	int score;
};

class Alignment{
private:
	string mstr1;
	string mstr2;
	int mmatch;
	int mmismatch;
	int mh;
	int mg;
public:	
	Alignment(int match=0, int mismatch=0, int h=0, int g=0);
	~Alignment();
	int readFasta(fstream fasta);
	int optimalGlobalAlignment();
	int optimalLocalAlignment();
};
