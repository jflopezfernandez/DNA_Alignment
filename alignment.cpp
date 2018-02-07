#include "alignment.h"

Alignment::Alignment(int match, int mismatch, int h, int g)
{
	mmatch = match;
	mmismatch = mismatch;
	mh = h;
	mg = g;
	mstr1 = "";
	mstr2 = "";
}
Alignment::Alignment(const Alignment &copy)
{
	mmatch = copy.mmatch;
	mmismatch = copy.mmismatch;
	mh = copy.mh;
	mg = copy.mg;
	mstr1 = "";
	mstr2 = "";
}
Alignment::~Alignment()
{

}
// file must be open 
int Alignment::readFasta(const char *fasta)
{
	int success = -1;
	char x = 0;
	ifstream fin;
	fin.open(fasta);
	cout << "readFasta()" << endl;
	
	if (fin.is_open())
	{
		success = 1;
//		cout << "reading " << endl;
		while(!fin.eof())
		{
			x = fin.get();
			cout << x;
		}
	}
	
	fin.close();
	return success;
}
int Alignment::optimalGlobalAlignment()
{
	return -1;
}
int Alignment::optimalLocalAlignment()
{
	return -1;
}
