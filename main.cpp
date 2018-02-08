#include "alignment.h"

void printAllArgs(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	int success = 0;
	ifstream fin;
	printAllArgs(argc, argv);
	Alignment align = Alignment();
	
	//<exe><inputfile w/ s1 s2><0:global, 1:local>[config file]
	if(argc == 3)
	{
		success == align.readFasta(argv[1]);
		success &= align.readConfig("parameters.config");
		align.optimalGlobalAlignment();
	}
	else if(argc > 3)
	{//config file specified 
		success == align.readFasta(argv[1]);
		success &= align.readConfig(argv[3]);
	}
	else
	{
		cout << "format: <exe name> <inputfile containing s1 s2> <0:global, 1:local> [parameter config file]" << endl;
	}
	cout << "success: ";
	success != 0? cout << "True" << endl: cout << "False" << endl;

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
