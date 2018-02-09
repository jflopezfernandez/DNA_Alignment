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
	cout << "s1 length: " << mstr1_len << endl << "s2 length: " << mstr2_len << endl;
	success &= setStringsFromFasta(fasta);
//	cout << "str1: \n" << mstr1 <<  str2: \n" << mstr2 << endl;
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
		for(int i = 0; fin.get(c) && i < mstr1_len; )
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
		for(int i = 0; fin.get(c) && i < mstr2_len; )
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
//	cout << "readConfig()" << endl;
	
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
		
		cout << "Scores: ";
		cout << "match=" << num_match << ", mismatch= " << num_mismatch;
		cout << ", h=" << num_h << ", g=" << num_g << " " << endl;

		mmatch = num_match;
		mmismatch = num_mismatch;
		mh = num_h;
		mg = num_g;
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
//	int aAlign[mstr1_len][mstr2_len] = { 0 };
	// use str lens to make the table 
	mMatrix = new DP_Cell*[mstr1_len];
	for (int i = 0; i < mstr1_len; i++)
		mMatrix[i] = new DP_Cell[mstr2_len];

	//init initial row/coloumn
	initFirstRowColumn();
	//go forward
	forwardComputationGlobal();
//	printMatrix();
	int optimal_score = maxOfThree(mMatrix[mstr1_len-1][mstr2_len-1].score_d,
					mMatrix[mstr1_len-1][mstr2_len-1].score_s,
					mMatrix[mstr1_len-1][mstr2_len-1].score_i);
//	cout << "optimal score" << optimal_score << endl;
	//retrace back 
	retraceGlobal();
	return optimal_score;
}
//aAlign[mstr1_len][mstr2_len]
int Alignment::initFirstRowColumn()
{
	//INT_MIN is min interger const from <climits>
	//init 0,0
	mMatrix[0][0].score_d = INT_MIN;
	mMatrix[0][0].score_s = 0;
	mMatrix[0][0].score_i = INT_MIN;
	//init i,0 where i > 0
	for(int i = 1; i < mstr1_len; i++)
	{
		mMatrix[i][0].score_d = mh + ( i * mg );
		mMatrix[i][0].score_s = INT_MIN;
		mMatrix[i][0].score_i = INT_MIN;
	}
	//init 0,j where j > 0
	for(int j = 1; j < mstr2_len; j++)
	{
		mMatrix[0][j].score_d = INT_MIN;
		mMatrix[0][j].score_s = INT_MIN;
		mMatrix[0][j].score_i = mh + ( j * mg );
	}
	return 0;
}
//row and columns must be initialized 
int Alignment::forwardComputationGlobal()
{
	int score_d_d, score_d_s, score_d_i; 
	int score_s_d, score_s_s, score_s_i; 
	int score_i_d, score_i_s, score_i_i; 

	//init just cuz 
	score_d_d = score_d_s = score_d_i = INT_MIN; 
	score_s_d = score_s_s = score_s_i = INT_MIN; 
	score_i_d = score_i_s = score_i_i = INT_MIN; 
	
	//start at 1, first row/column already done
	for(int i = 1; i < mstr1_len; i++)
	{
		for(int j = 1; j < mstr2_len; j++)
		{
			// calculate d score 
			if(mMatrix[i-1][j].score_d != INT_MIN) score_d_d = mMatrix[i-1][j].score_d + mg;
			if(mMatrix[i-1][j].score_s != INT_MIN) score_d_s = mMatrix[i-1][j].score_s + mg + mh; // start of gap
			if(mMatrix[i-1][j].score_i != INT_MIN) score_d_i = mMatrix[i-1][j].score_i + mg + mh; // start of gap
			// calculate s score 
			if(mMatrix[i-1][j-1].score_d != INT_MIN) score_s_d = mMatrix[i-1][j-1].score_d + S(i,j);
			if(mMatrix[i-1][j-1].score_s != INT_MIN) score_s_s = mMatrix[i-1][j-1].score_s + S(i,j);
			if(mMatrix[i-1][j-1].score_i != INT_MIN) score_s_i = mMatrix[i-1][j-1].score_i + S(i,j);
			// calculate i score 
			if(mMatrix[i][j-1].score_d != INT_MIN) score_i_d = mMatrix[i][j-1].score_d + mg + mh; // start of gap
			if(mMatrix[i][j-1].score_s != INT_MIN) score_i_s = mMatrix[i][j-1].score_s + mg + mh; // start of gap
			if(mMatrix[i][j-1].score_i != INT_MIN) score_i_i = mMatrix[i][j-1].score_i + mg;
			// find maxes
			mMatrix[i][j].score_d = maxOfThree(score_d_d, score_d_s, score_d_i);
			mMatrix[i][j].score_s = maxOfThree(score_s_d, score_s_s, score_s_i);
			mMatrix[i][j].score_i = maxOfThree(score_i_d, score_i_s, score_i_i);
		}
	}
	return 1;
}
int Alignment::maxOfThree(int n1, int n2, int n3)
{
	int max = n1;
	if(max < n2)
		max = n2;
	if(max < n3)
		max = n3;
	return max;
}
//returns match/mismatch score of a cell 
int Alignment::S(int i, int j)
{
	int score = 0;
	if( mstr1[i] == mstr2[j] )
	{
		score = mmatch;
	}
	else
	{
		score = mmismatch;
	}
	return score;
}
//str(1/2)_len and matrix must be intialized 
int Alignment::printMatrix()
{
	for(int i = 0; i < mstr1_len; i++)
	{
		for(int j = 0; j < mstr2_len; j++)
		{
			cout << "[" << mMatrix[i][j].score_d;
			cout << "," << mMatrix[i][j].score_s;
			cout << "," << mMatrix[i][j].score_i << "] ";
		}
		cout << endl;
	}
}
// things must be initialized 
// THis will print to stdout 
void Alignment::retraceGlobal()
{
	//start at last cell and move to the cell with the max score 
	int i = mstr1_len - 1;
	int j = mstr2_len - 1;
	int cell_max = 0;
	char aStr1[mstr1_len + mstr2_len] = { 0 };
	char aStr2[mstr1_len + mstr2_len] = { 0 };
	int aStr_index = 0;
	int aStr2_index = 0;
	int matches = 0, mismatches = 0, gaps = 0, opening_gaps = 0;
	bool bGap = false;

	
	// record the aligned strings 
	while ( i != -1 || j != -1 )
	{
		//find whether D,S, or I was the max 
		cell_max = maxOfThree(mMatrix[i][j].score_d,
					 mMatrix[i][j].score_s,
					 mMatrix[i][j].score_i); 
		if(cell_max == mMatrix[i][j].score_d )
		{//deletion
			aStr1[aStr_index] = mstr1[i];
			aStr2[aStr_index] = '-';
			i--;
			if (bGap == false)
				opening_gaps++;
			bGap = true;
			gaps++;
		}
		else if( cell_max == mMatrix[i][j].score_s)
		{//subsitution 	
			aStr1[aStr_index] = mstr1[i];
			aStr2[aStr_index] = mstr2[j];
			i--;
			j--;
			bGap = false;
			if (aStr1[aStr_index] == aStr2[aStr_index])
				matches++;
			else
				mismatches++;
		}
		else 
		{//insertion 	
			aStr1[aStr_index] = '-';
			aStr2[aStr_index] = mstr2[j];
			j--;
			if (bGap == false)
				opening_gaps++;
			bGap = true;
			gaps++;
		}

		aStr_index++;

		if ( i < -1 || j < -1)
		{
			cout << "i or j is < 0... oops... abort..." << endl;
			return;
		}
	}
	aStr_index--; //I want it to be on the last character 
	//print that alignment!
	int count = 1, count2 = 1;
	bool extra_tab = false;
	while(aStr_index >= 0)
	{
		//print 60 or more charaters of str1
		if((count/10000) >= 1)
			extra_tab = true;
		cout << "s1 " << count << '\t';
		for(int num = 0; num < 60 && (aStr_index - num) >= 0; num++)
		{
			cout << aStr1[aStr_index - num];
			if(aStr1[aStr_index - num] != '-') 
				count++;
		}
		cout << ' ' << (count - 1) << endl;
		
		//print  the | and ' '
		cout << '\t';
		if ( extra_tab )
			cout << '\t';
		for(int num = 0; num < 60 && (aStr_index - num) >= 0; num++)
		{
			if(aStr1[aStr_index - num] == aStr2[aStr_index - num])
				cout << '|';
			else
				cout << ' ';
		}
		cout << endl;
		

		//print 60 or more charaters of str2
		cout << "s2 " << count2 << '\t';
		if ( (count2 /10000) < 1 && extra_tab )
			cout << '\t';
		for(int num = 0; num < 60 && (aStr_index - num) >= 0; num++)
		{
			cout << aStr2[aStr_index - num];
			if(aStr2[aStr_index - num] != '-')
				count2++;
		}
		cout << ' ' << (count2 - 1) << endl << endl;
		
		aStr_index -= 60;
	}	
	//print them stats 
	int len = 0;
	mstr1_len > mstr2_len? len = mstr1_len: len = mstr2_len;
	cell_max = maxOfThree(mMatrix[mstr1_len-1][mstr2_len-1].score_d,
				 mMatrix[mstr1_len-1][mstr2_len-1].score_s,
				 mMatrix[mstr1_len-1][mstr2_len-1].score_i); 
	cout << "Global alignment score: " << cell_max <<  endl;
	cout << "Number of: matches=" << matches << ", mismatches=" << mismatches;
	cout << ", gaps=" << gaps << ", opening gaps=" << opening_gaps << endl;
	cout << "Identities= " << matches << '/' << len << " (";
	cout << ((matches*100)/len) << "%), Gaps =" << gaps << '/' << len;
	cout << " (" << ((gaps*100)/len) << "%)" << endl;
}
int Alignment::optimalLocalAlignment()
{
	// use str lens to make the table

	// As a note... Both alignment functions 
	// allocate a lot of memory and since I am 
	// only running this one or the other one 
	// once per execution, I don't actually 
	// bother to deallocate any of it. 
	// So, beware. 
	mMatrix = new DP_Cell*[mstr1_len];
	for (int i = 0; i < mstr1_len; i++)
		mMatrix[i] = new DP_Cell[mstr2_len];

	//init initial row/coloumn
	initFirstRowColumnLocal();

	//go forward
	forwardComputationLocal();
//	printMatrix();

	int optimal_score = maxOfThree(mMatrix[mstr1_len-1][mstr2_len-1].score_d,
					mMatrix[mstr1_len-1][mstr2_len-1].score_s,
					mMatrix[mstr1_len-1][mstr2_len-1].score_i);
//	cout << "optimal score" << optimal_score << endl;

	//retrace back 
	retraceLocal();
	return optimal_score;
}
int Alignment::initFirstRowColumnLocal()
{
	//INT_MIN is min interger const from <climits>
	//init i,0 where i >= 0
	for(int i = 0; i < mstr1_len; i++)
	{
		mMatrix[i][0].score_d = 0;
		mMatrix[i][0].score_s = 0;
		mMatrix[i][0].score_i = 0;
	}
	//init 0,j where j >= 0
	for(int j = 0; j < mstr2_len; j++)
	{
		mMatrix[0][j].score_d = 0;
		mMatrix[0][j].score_s = 0;
		mMatrix[0][j].score_i = 0;
	}
	return 0;
}
//row and columns must be initialized 
int Alignment::forwardComputationLocal()
{
	int score_d_d, score_d_s, score_d_i; 
	int score_s_d, score_s_s, score_s_i; 
	int score_i_d, score_i_s, score_i_i; 
	int cell_max = 0;
	mlocal_max = 0;  //member variable
	mlocal_max_i = 0;
	mlocal_max_j = 0;

	//init just cuz 
	score_d_d = score_d_s = score_d_i = 0; 
	score_s_d = score_s_s = score_s_i = 0; 
	score_i_d = score_i_s = score_i_i = 0; 
	
	//start at 1, first row/column already done
	for(int i = 1; i < mstr1_len; i++)
	{
		for(int j = 1; j < mstr2_len; j++)
		{
			// calculate d score 
			if(mMatrix[i-1][j].score_d != INT_MIN) score_d_d = mMatrix[i-1][j].score_d + mg;
			if(mMatrix[i-1][j].score_s != INT_MIN) score_d_s = mMatrix[i-1][j].score_s + mg + mh; // start of gap
			if(mMatrix[i-1][j].score_i != INT_MIN) score_d_i = mMatrix[i-1][j].score_i + mg + mh; // start of gap
			// calculate s score 
			if(mMatrix[i-1][j-1].score_d != INT_MIN) score_s_d = mMatrix[i-1][j-1].score_d + S(i,j);
			if(mMatrix[i-1][j-1].score_s != INT_MIN) score_s_s = mMatrix[i-1][j-1].score_s + S(i,j);
			if(mMatrix[i-1][j-1].score_i != INT_MIN) score_s_i = mMatrix[i-1][j-1].score_i + S(i,j);
			// calculate i score 
			if(mMatrix[i][j-1].score_d != INT_MIN) score_i_d = mMatrix[i][j-1].score_d + mg + mh; // start of gap
			if(mMatrix[i][j-1].score_s != INT_MIN) score_i_s = mMatrix[i][j-1].score_s + mg + mh; // start of gap
			if(mMatrix[i][j-1].score_i != INT_MIN) score_i_i = mMatrix[i][j-1].score_i + mg;
			// find maxes
			mMatrix[i][j].score_d = maxOfThree(score_d_d, score_d_s, score_d_i);
			mMatrix[i][j].score_s = maxOfThree(score_s_d, score_s_s, score_s_i);
			mMatrix[i][j].score_i = maxOfThree(score_i_d, score_i_s, score_i_i);
			// zero out any negative scores 
			if(mMatrix[i][j].score_d < 0) mMatrix[i][j].score_d = 0;
			if(mMatrix[i][j].score_s < 0) mMatrix[i][j].score_s = 0;
			if(mMatrix[i][j].score_i < 0) mMatrix[i][j].score_i = 0;
			
			// also I want to know the max cell 
			cell_max = maxOfThree(mMatrix[i][j].score_d, mMatrix[i][j].score_s, mMatrix[i][j].score_i);
			if (cell_max > mlocal_max)
			{
				mlocal_max = cell_max;
				mlocal_max_i = i;
				mlocal_max_j = j;
			}
		}
	}
	return 1;
}
void Alignment::retraceLocal()
{
	//start at max cell
	int i = mlocal_max_i;
	int j = mlocal_max_j;
	int cell_max = 0;
	char aStr1[mstr1_len + mstr2_len] = { 0 };
	char aStr2[mstr1_len + mstr2_len] = { 0 };
	int aStr_index = 0;
	int aStr2_index = 0;
	int matches = 0, mismatches = 0, gaps = 0, opening_gaps = 0;
	bool bGap = false;

	// record the aligned strings 
	while ( i != -1 || j != -1 )
	{
		//find whether D,S, or I was the max 
		cell_max = maxOfThree(mMatrix[i][j].score_d,
					 mMatrix[i][j].score_s,
					 mMatrix[i][j].score_i); 
		if(cell_max == 0)
		{
			break;
		}
		if(cell_max == mMatrix[i][j].score_d )
		{//deletion
			aStr1[aStr_index] = mstr1[i];
			aStr2[aStr_index] = '-';
			i--;
			if (bGap == false)
				opening_gaps++;
			bGap = true;
			gaps++;
		}
		else if( cell_max == mMatrix[i][j].score_s)
		{//subsitution 	
			aStr1[aStr_index] = mstr1[i];
			aStr2[aStr_index] = mstr2[j];
			i--;
			j--;
			bGap = false;
			if (aStr1[aStr_index] == aStr2[aStr_index])
				matches++;
			else
				mismatches++;
		}
		else 
		{//insertion 	
			aStr1[aStr_index] = '-';
			aStr2[aStr_index] = mstr2[j];
			j--;
			if (bGap == false)
				opening_gaps++;
			bGap = true;
			gaps++;
		}

		aStr_index++;

		if ( i < -1 || j < -1)
		{
			cout << "i or j is < 0... oops... abort..." << endl;
			return;
		}
	}
	aStr_index--; //I want it to be on the last character 
	//print that alignment!
	int count = 1, count2 = 1;
	bool extra_tab = false;
	while(aStr_index >= 0)
	{
		if ( (count / 10000) >= 1)
			extra_tab = true;
		//print 60 or more charaters of str1
		cout << "s1 " << count << '\t';
		for(int num = 0; num < 60 && (aStr_index - num) >= 0; num++)
		{
			cout << aStr1[aStr_index - num];
			if(aStr1[aStr_index - num] != '-') 
				count++;
		}
		cout << ' ' << (count - 1) << endl;
		
		//print  the | and ' '
		cout << '\t';
		if ( extra_tab )
			cout << '\t';
		for(int num = 0; num < 60 && (aStr_index - num) >= 0; num++)
		{
			if(aStr1[aStr_index - num] == aStr2[aStr_index - num])
				cout << '|';
			else
				cout << ' ';
		}
		cout << endl;
		

		//print 60 or more charaters of str2
		cout << "s2 " << count2 << '\t';
		if ( (count2 / 10000) < 1 && extra_tab)
			cout << '\t';
		for(int num = 0; num < 60 && (aStr_index - num) >= 0; num++)
		{
			cout << aStr2[aStr_index - num];
			if(aStr2[aStr_index - num] != '-')
				count2++;
		}
		cout << ' ' << (count2 - 1) << endl << endl;
		
		aStr_index -= 60;
	}	
	//print them stats 
	int len = 0;
	count > count2? len = (count - 1): len = (count2 - 1);

	cout << "Local alignment score: " << mlocal_max <<  endl;
	cout << "Number of: matches=" << matches << ", mismatches=" << mismatches;
	cout << ", gaps=" << gaps << ", opening gaps=" << opening_gaps << endl;
	cout << "Identities= " << matches << '/' << len << " (";
	cout << ((matches*100)/len) << "%), Gaps =" << gaps << '/' << len;
	cout << " (" << ((gaps*100)/len) << "%)" << endl;

}
