#include "alignment.h"

Alignment::Alignment(int match, int mismatch, int h, int g)
{
	mmatch = match;
	mmismatch = mismatch;
	mh = h;
	mg = g;
	mstr1 = "";
	mstr2 = "";
	mstr1_len = 0;
	mstr2_len = 0;
}
Alignment::Alignment(const Alignment &copy)
{
	mmatch = copy.mmatch;
	mmismatch = copy.mmismatch;
	mh = copy.mh;
	mg = copy.mg;
	mstr1 = "";
	mstr2 = "";
	mstr1_len = 0;
	mstr2_len = 0;
}
Alignment::~Alignment()
{

} 
int Alignment::readFasta(const char *fasta)
{
	int success = 0;
	success  = setLengthsFromFasta(fasta);
	cout << "len1: " << mstr1_len << ", len2: " << mstr2_len << endl;
	success &= setStringsFromFasta(fasta);
	cout << "str1: \n" << mstr1 << ", str2: \n" << mstr2 << endl;
	return success;
}
int Alignment::setLengthsFromFasta(const char *fasta)
{	
	int success = 0;
	ifstream fin;
	char line[1024];
	char c = 0;

	mstr1_len = 0;
	mstr2_len = 0;

	fin.open(fasta);
	if (fin.is_open())
	{
		//ignore the header
		fin.ignore(10, '>');
		fin.getline(line, 1024);
		//count str1_len
		while(!fin.eof() && fin.get(c))
		{
			if(c == '>')
				break;
			if(c != '\n')
				mstr1_len++;
		}
		//ignore the header 
		fin.getline(line, 1024);
		//count str2_len
		while(!fin.eof() && fin.get(c))
		{
			if(c != '\n')
				mstr2_len++;
		}
	}
	else
	{
		cout << "file " << fasta << "could not open..." << endl;
	}
	fin.close();
	return success;
}
// must have correct file format and mstr(1/2)_len must be initialized
// by setLengthsFromFasta()
int Alignment::setStringsFromFasta(const char *fasta)
{	
	int success = 0;
	ifstream fin;
	char line[1024];
	char c = 0;

	mstr1.resize(mstr1_len + 1);
	mstr2.resize(mstr2_len + 1);

	fin.open(fasta);
	if (fin.is_open())
	{
		//ignore the header
		fin.ignore(10, '>');
		fin.getline(line, 1024);
		//count str1_len
		for(int i = 0;/* !fin.eof() &&*/ fin.get(c) && i < mstr1_len; )
		{
			if(c == '>')
				break;
			if(c != '\n')
			{
				mstr1[i] = c;
				i++;
			}
		}
		mstr1[mstr1_len + 1] = '\0';
		//ignore the header 
		fin.getline(line, 1024);
		fin.getline(line, 1024);
//		cout << "line=" << line << endl;
		//count str2_len
		for(int i = 0;/* !fin.eof() &&*/ fin.get(c) && i < mstr2_len; )
		{
			if(c != '\n')
			{
				mstr2[i] = c;
				i++;
			}
		}
		mstr2[mstr2_len + 1] = '\0';
	}
	else
	{
		cout << "file " << fasta << "could not open..." << endl;
	}
	fin.close();
	return success;
}
int Alignment::readConfig(const char *config)
{
	int success = -1;
	char x = 0;
	int num_match = 0;
	int num_mismatch = 0;
	int num_h = 0;
	int num_g = 0;
	ifstream fin;
	fin.open(config);
	cout << "readConfig()" << endl;
	
	if (fin.is_open())
	{
		success = 1;
		fin.ignore(20, ' ');
		fin >> num_match;
		fin.ignore(20, ' ');
		fin >> num_mismatch;
		fin.ignore(20, ' ');
		fin >> num_h;
		fin.ignore(20, ' ');
		fin >> num_g;
		
		cout << num_match << " " << num_mismatch << " ";
		cout << num_h << " " << num_g << " " << endl;
	}
	else
	{
		cout << config << " read failed..." << endl;
	}
	
	fin.close();
	return success;
}
int Alignment::optimalGlobalAlignment()
{
	int aAlign[mstr1_len][mstr2_len] = { 0 };
	// use str lens to make the table 
	mMatrix = new DP_Cell*[mstr1_len];
	for (int i = 0; i < mstr1_len; i++)
		mMatrix[i] = new DP_Cell[mstr2_len];

	//init initial row/coloumn
	initFirstRowColumn();
	//go forward 
	//retrace back 
	return -1;
}
//aAlign[mstr1_len][mstr2_len]
int Alignment::initFirstRowColumn()
{
	//INT_MIN is min interger const from <climits>
	//init 0,0
	mMatrix[0][0].score_d = 0;
	mMatrix[0][0].score_s = 0;
	mMatrix[0][0].score_i = 0;
	//init i,0
	for(int i = 1; i < mstr1_len; i++)
	{
		mMatrix[i][0].score_d = mh + ( i * mg );
		mMatrix[i][0].score_s = INT_MIN;
		mMatrix[i][0].score_i = INT_MIN;
	}
	//init 0,j
	for(int j = 0; j < mstr2_len; j++)
	{
		mMatrix[0][j].score_d = INT_MIN;
		mMatrix[0][j].score_s = INT_MIN;
		mMatrix[0][j].score_i = mh + ( j * mg );
	}
	return 0;
}
int Alignment::optimalLocalAlignment()
{
	return -1;
}
