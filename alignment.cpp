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
		while(!fin.eof())
		{
			//read until > get the header line 
			//read s1 
			//read until > get the header line
			//read s2 
			x = fin.get();
			cout << x;
		}
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
	return -1;
}
int Alignment::optimalLocalAlignment()
{
	return -1;
}
