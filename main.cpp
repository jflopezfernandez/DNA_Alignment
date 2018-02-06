#include "alignment.h"

void printAllArgs(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	printAllArgs(argc, argv);
	Alignment align();	
}

void printAllArgs(int argc, char *argv[])
{
	for(int i = 0; i < argc; i++)
	{
		cout << argv[i] << " ";
	}
	cout << endl;
}
