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
int Alignment::readFasta(const char *fasta)
{
	int success = -1;
	int str_start = 0;
	int str_len = 0;
	char c = 0;
	char line[1024];
	ifstream fin;
	fin.open(fasta);
	cout << "readFasta()" << endl;
	
	if (fin.is_open())
	{
		success = 1;
		cout << "It's opened!" << endl;
		while(!fin.eof())
		{
			fin.ignore(10, '>');
			cout << "did the ignore" << endl;
			fin.getline(line, 1024);
			cout << "did the getline" << endl;
			cout << line << endl;
			
			str_start = fin.tellg();
			cout << "start of string: " << str_start << endl;
		
			while( !fin.eof() && fin.get(c))
			{
				if(c == '\n')
				{
					if(!fin.eof() && fin.get(c))
					{
						if (c == '\n')
							break;
						else
							str_len++;
					}
				}
				else
				{
					str_len++;
				}
			}
			if (c != '\n')
				return 0;
			cout << "length of string: " << str_len << endl;
			
//			fin.getline(line, 512);
//			cout << line << endl;

			fin.seekg(str_start, fin.beg);
			mstr1.resize(str_len + 1);
			for (int i = 0; i < str_len; i++)
			{
				fin.get(c);
				cout << c;
				if ( c != '\n' )
					mstr1[i] = c;
				else
					i--;
			}
			mstr1[str_len + 1] = '\0';
			cout << endl;


			break;	
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
