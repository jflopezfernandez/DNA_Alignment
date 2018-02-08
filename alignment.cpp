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
/*	int str_start = 0;
//	int str_len = 0;
	char c = 0;
	char line[1024];
	ifstream fin;
	fin.open(fasta);
	cout << "readFasta()" << endl;
	
	if (fin.is_open())
	{
		success = 1;
		cout << "It's opened!" << endl;
	
		fin.ignore(10, '>');
//		cout << "did the ignore" << endl;
		fin.getline(line, 1024);
//		cout << "did the getline" << endl;
//		cout << line << endl;
			
		str_start = fin.tellg();
		cout << "start of string: " << str_start << endl;
		
		//length of str1 
		while( !fin.eof() && fin.get(c))
		{
			if(c == '\n')
			{
				if(!fin.eof() && fin.get(c))
				{
					if (c == '\n')
						break;
					else
						mstr1_len++;
				}
			}
			else
			{
				mstr1_len++;
			}
		}
		cout << "length of string 1: " << mstr1_len << endl;
		
//		fin.getline(line, 512);
//		cout << line << endl;
		
		//copy the seq over to the str skipping \n
		fin.seekg(str_start, fin.beg);
		mstr1.resize(mstr1_len + 1);
		for (int i = 0; i < mstr1_len; i++)
		{
			fin.get(c);
//			cout << c;
			if ( c != '\n' )
				mstr1[i] = c;
			else
				i--;
		}
		mstr1[mstr1_len + 1] = '\0';
		cout << mstr1 << endl;
		
		fin.ignore(3, '>');
		fin.getline(line, 1024);
		cout << endl << line << endl;

		// get str2 length 
		while(fin.get(c))
		{
			cout << c << endl;
			if(c != '\n')
				break;
		}
//		fin.get(c);
		c == '\n'? cout << "Fuck that" << endl: cout << "Shit... now what" << endl;
		str_start = fin.tellg();
		str_start--;
		cout << "start of 2nd string: " << str_start << endl;
		while(!fin.eof() && fin.get(c))
		{
			if(c != '\n')
				mstr2_len++;
		}

		mstr2.resize(mstr2_len + 1);
		cout << "length of 2nd string: " << mstr2_len << endl;
		
		//copy seq to str2
		
		fin.seekg(str_start, fin.beg);
//		fin.getline(line, 1024);
//		cout << "line at start: " << line << endl;
		for (int i = 0; i < mstr2_len; i++)
		{
			fin.get(c);
			cout << c;
			if ( c != '\n' )
				mstr2[i] = c;
			else
				i--;
		}
		mstr2[mstr2_len + 1] = '\0';
		cout << mstr2 << endl;
		
	}
	
	fin.close();*/
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
		for(int i = 0; !fin.eof() && fin.get(c) && i < mstr1_len; )
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
		cout << "line=" << line << endl;
		//count str2_len
		for(int i = 0; !fin.eof() && fin.get(c) && i < mstr2_len; )
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
	return -1;
}
int Alignment::optimalLocalAlignment()
{
	return -1;
}
