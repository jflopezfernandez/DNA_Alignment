/*
* Michelle Farr
*/
#include <iostream>
#include <string>
#include <fstream>
#include <climits>

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::ifstream;

struct DP_Cell{
	int score_s;	//   
	int score_d;	//    S D
	int score_i;	//    I[ ]
};

class Alignment{
private:
	string mstr1;
	string mstr2;
	int mstr1_len;
	int mstr2_len;
	int mmatch;
	int mmismatch;
	int mh;
	int mg;
	DP_Cell **mMatrix;
	
	int setLengthsFromFasta(const char *fasta);
	int setStringsFromFasta(const char *fasta);
	int initFirstRowColumn();
	int forwardComputationGlobal();
	int maxOfThree(int n1, int n2, int n3);
	int S(int i, int j);
	int printMatrix();
	void retraceGlobal();
public:	
	Alignment(int match=0, int mismatch=0, int h=0, int g=0);
	Alignment(const Alignment &copy);
	~Alignment();
	int readFasta(const char *fasta);
	int readConfig(const char *config);
	int optimalGlobalAlignment();
	int optimalLocalAlignment();
};
