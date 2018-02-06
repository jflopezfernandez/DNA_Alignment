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
Alignment::~Alignment()
{

}
int Alignment::readFasta(fstream fasta)
{
	return -1;
}
int Alignment::optimalGlobalAlignment()
{
	return -1;
}
int Alignment::optimalLocalAlignment()
{
	return -1;
}
