#include "alignment.h"

void printAllArgs(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	int success = 0;
	ifstream fin;
	printAllArgs(argc, argv);
	cout <<  endl;
	Alignment align = Alignment();
	
	//<exe><inputfile w/ s1 s2><0:global, 1:local>[config file]
	if(argc == 3)
	{
		success == align.readFasta(argv[1]);
		success &= align.readConfig("parameters.config");
		if(argv[2][0] == '0')
			align.optimalGlobalAlignment();
		else 
			align.optimalLocalAlignment();
	}
	else if(argc > 3)
	{//config file specified 
		success == align.readFasta(argv[1]);
		success &= align.readConfig(argv[3]);
		if(argv[2][0] == '0')
			align.optimalGlobalAlignment();
		else
			align.optimalLocalAlignment();
	}
	else
	{
		cout << "format: <exe name> <inputfile containing s1 s2> <0:global, 1:local> [parameter config file]" << endl;
	}
	return 0;
}

void printAllArgs(int argc, char *argv[])
{
	for(int i = 0; i < argc; i++)
	{
		cout << argv[i] << " ";
	}
	cout << endl;
}
